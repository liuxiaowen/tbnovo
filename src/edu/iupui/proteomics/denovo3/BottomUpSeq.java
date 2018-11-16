package edu.iupui.proteomics.denovo3;

import java.util.ArrayList;
import java.util.Arrays;

import edu.iupui.proteomics.base.seq.ResSeq;
import edu.iupui.proteomics.base.theosp.BpSpec;
import edu.iupui.proteomics.spec.peak.PeakTolerance;
import edu.iupui.proteomics.spec.prmsp.PrmPeak;
import edu.iupui.proteomics.spec.sp.Ms;

public class BottomUpSeq {
	// // bottom-up de novo resulting sequences
	public BpSpec seq;
	public double nGap;
	public double cGap;
	// bottom-up prefix residue mass spectra
	private double seqPrms[];
	// best shifts and scores for bottom-up mass spectra
	private double bestShift = -1;
	private double bestScore = -1;
	private double weight = 1;
	private int matchStart = -1;
	private int matchEnd = -1;
	
	/* ms/ms information */
	private String fileName;
	private int scan = -1;
	private double denovoScore;
	private Ms<PrmPeak> prmSp;
	
	public BottomUpSeq(BpSpec seq, double nGap, double cGap) {
		this.seq = seq;
		seqPrms = seq.getExtBMasses();
		this.nGap = nGap;
		this.cGap = cGap;
		double basePrms[] = seq.getExtBMasses();
		ArrayList<Double> prmList = new ArrayList<Double>();
		prmList.add(0.0);
		if (nGap != 0.0) {
			prmList.add(nGap);
		}
		for (int i = 1; i < basePrms.length; i++) {
			prmList.add(nGap + basePrms[i]);
		}
		if (cGap != 0.0) {
			prmList.add(nGap+basePrms[basePrms.length-1] + cGap);
		}
		seqPrms = new double[prmList.size()];
		for (int i = 0; i < prmList.size(); i++) {
			seqPrms[i] = prmList.get(i);
		}
	}
	
	public double getNGap() {
		return nGap;
	}

	public BottomUpSeq(BpSpec seq, double nGap, double cGap, String fileName, int scan, double denovoScore) {
		this.seq = seq;
		this.scan = scan;
		this.fileName = fileName;
		this.setDenovoScore(denovoScore);
		this.nGap = nGap;
		this.cGap = cGap;
		double basePrms[] = seq.getExtBMasses();
		ArrayList<Double> prmList = new ArrayList<Double>();
		prmList.add(0.0);
		if (nGap != 0.0) {
			prmList.add(nGap);
		}
		for (int i = 1; i < basePrms.length; i++) {
			prmList.add(nGap + basePrms[i]);
		}
		if (cGap != 0.0) {
			prmList.add(nGap+basePrms[basePrms.length-1] + cGap);
		}
		seqPrms = new double[prmList.size()];
		for (int i = 0; i < prmList.size(); i++) {
			seqPrms[i] = prmList.get(i);
		}
	}
	
	private void updateStartEnd (double masses[], PeakTolerance peakTolerance) {
		matchStart = -1;
		matchEnd = -1;
		double shiftPrms[] = getShiftedPrms();
		for (int i = 0; i < shiftPrms.length; i++) {
			boolean found = false;
			for (int j = 0; j < masses.length; j++) {
				double diff = Math.abs(shiftPrms[i] - masses[j]);
				if (diff <= peakTolerance.compStrictErrorTole(masses[j])) {
					found = true;
				}
			}
			if (found) {
				matchStart = i;
				break;
			}
		}
		for (int i = shiftPrms.length - 1; i >= 0; i--) {
			boolean found = false;
			for (int j = 0; j < masses.length; j++) {
				double diff = Math.abs(shiftPrms[i] - masses[j]);
				if (diff <= peakTolerance.compStrictErrorTole(masses[j])) {
					found = true;
				}
			}
			if (found) {
				matchEnd = i;
				break;
			}
		}
	}

	public void findBestShift(double masses[], PeakTolerance peakTolerance) {
		bestShift = -1;
		bestScore = -1;
		for (int i = 0; i < masses.length; i++) {
			for (int j = 0; j < seqPrms.length; j++) {
				double shift = masses[i] - seqPrms[j];
				double score = compDiagScr(masses, seqPrms, shift,
						peakTolerance);
				if (score > bestScore) {
					bestScore = score;
					bestShift = shift;
				}
			}
		}
		updateStartEnd(masses, peakTolerance);
		return;
	}
	
	public void findLimitedBestShift(double masses[], double shifts[],
			PeakTolerance peakTolerance) {
		bestShift = -1;
		bestScore = -1;
		for (int i = 0; i < shifts.length; i++) {
			/* WORKING ON THIS */
			double score = compDiagScr(masses, seqPrms, shifts[i],
					peakTolerance);
//			score = 0.5 * compDiagScr(masses, seqPrms, shifts[i]+MassConstant.getIsotopeMass(), peakTolerance);
//			score = 0.5 * compDiagScr(masses, seqPrms, shifts[i]-MassConstant.getIsotopeMass(), peakTolerance);
			if (score > bestScore) {
				bestScore = score;
				bestShift = shifts[i];
			}
		}
	}

	public void findEndBestShift(double masses[], PeakTolerance peakTolerance) {
		bestShift = -1;
		bestScore = -1;
		for (int i = 0; i < masses.length; i++) {
			double shift = masses[i];
			double score = compDiagScr(masses, seqPrms, shift, peakTolerance);
			if (score > bestScore) {
				bestScore = score;
				bestShift = shift;
			}
			shift = masses[i] - seqPrms[seqPrms.length-1];
			score = compDiagScr(masses, seqPrms, shift, peakTolerance);
			if (score > bestScore) {
				bestScore = score;
				bestShift = shift;
			}
		}
		updateStartEnd(masses, peakTolerance);
		return;
	}

	protected static double compDiagScr(double masses[], double prms[],
			double center, PeakTolerance peakTolerance) {
		int i = 0;
		int j = 0;
		double scr = 0;
		while (i < masses.length && j < prms.length) {
			double distance = masses[i] - prms[j];
			if (Math.abs(center - distance) <= peakTolerance
					.compStrictErrorTole(masses[i])) {
				scr += 1.0;
				i++;
				j++;
			}
			/*
			 * we use 1 here since the difference between consecutive mass_b is
			 * at least 50. and some mass_a[i+1] may have large error tolerance
			 * than mass_a[i]
			 */
			if (distance > center + 1) {
				j++;
			} else {
				i++;
			}
		}
		return scr;
	}

	public double getBestScore() {
		return bestScore;
	}

	public double getNormalizedScore() {
		return bestScore / seqPrms.length;
	}

	public double getBestShift() {
		return bestShift;
	}

	public void setBestShift(double bestShift) {
		this.bestShift = bestShift;
	}

	public double getShiftedPrecMass() {
		return bestShift + seqPrms[seqPrms.length - 1];
	}

	public double[] getShiftedPrms() {
		double shiftedPrms[] = new double[seqPrms.length];
		for (int i = 0; i < seqPrms.length; i++) {
			shiftedPrms[i] = seqPrms[i] + bestShift;
		}
		return shiftedPrms;
	}
	
	public double[] getCysShiftedPrms(int pos) {
		double shiftedPrms[] = new double[seqPrms.length];
		for (int i = 0; i < seqPrms.length; i++) {
			if (i < pos) {
				shiftedPrms[i] = seqPrms[i] + bestShift;
			}
			else {
				shiftedPrms[i] = seqPrms[i] + bestShift + 57.021464;
			}
		}
		return shiftedPrms;
	}

	public double[] getShiftedTagPrms() {
		if (matchStart < 0) {
			return new double[0];
		}
		double shiftedPrms[] = new double[matchEnd - matchStart + 1];
		for (int i = matchStart; i <= matchEnd; i++) {
			shiftedPrms[i - matchStart] = seqPrms[i] + bestShift;
		}
		return shiftedPrms;
	}

	public BpSpec getBpSpec() {
		return seq;
	}

	public void output() {
		System.out.println(getBpSpec().getResSeq().getAcidString() + " shift "
				+ getBestShift() + " score " + getBestScore() + " percentage "
				+ ((double) (getBestScore()) / getShiftedPrms().length));
	}

	public double getWeight() {
		return weight;
	}

	public void setWeight(double weight) {
		this.weight = weight;
	}

	public int getScan() {
		return scan;
	}

	public String getFileName() {
		return fileName;
	}

	public void setFileName(String fileName) {
		this.fileName = fileName;
	}

	public Ms<PrmPeak> getPrmSp() {
		return prmSp;
	}

	public void setPrmSp(Ms<PrmPeak> prmSp) {
		this.prmSp = prmSp;
	}

	public double getDenovoScore() {
		return denovoScore;
	}

	public void setDenovoScore(double denovoScore) {
		this.denovoScore = denovoScore;
	}
	
	public int getCysteineNum() {
		ResSeq resSeq = seq.getResSeq();
		int n = 0;
		for (int i = 0; i < resSeq.getLen(); i++) {
			if (resSeq.getRes(i).getAcid().getOneLetter().equals("C")) {
				n++; 
			}
		}
		return n;
	}
	
	public double getSimilarResSumMass(double theoMass) {
		int n = getCysteineNum();
		double bestDiff = Double.MAX_VALUE;
		double bestMass = 0;
		double pepMass = getBpSpec().getResSeq().getResMassSum();
		for (int i = 0; i <= n; i++) {
			double m = pepMass + i * 57.021464;
			double d = Math.abs(m - theoMass);
			if (d < bestDiff) {
				bestDiff = d;
				bestMass = m;
			}
		}
		return bestMass;
	}
	
	public  static BottomUpSeq[] removeSimilarMass(BottomUpSeq seqs[]) {
		ArrayList<BottomUpSeq> shortList = new ArrayList<BottomUpSeq>();
		Arrays.sort(seqs, new DenovoScoreComparator());
		for (int i = 0; i < seqs.length; i++) {
			boolean found = false;
			for (int j = 0; j < shortList.size(); j++) {
				double mass1 = seqs[i].getBpSpec().getResSeq().getResMassSum();
				double mass2 = shortList.get(j).getBpSpec().getResSeq().getResMassSum();
				if (Math.abs(mass1-mass2) <= 0.05) {
					found = true;
				}
			}
			if (!found) {
				shortList.add(seqs[i]);
			}
		}
		return shortList.toArray(new BottomUpSeq[0]);
	}
}
