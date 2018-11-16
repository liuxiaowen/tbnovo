package edu.iupui.proteomics.denovo3;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

import org.apache.log4j.Logger;

import edu.iupui.proteomics.base.residue.MassConstant;
import edu.iupui.proteomics.base.residue.Res;
import edu.iupui.proteomics.base.residue.ResList;
import edu.iupui.proteomics.base.residue.ResListFactory;
import edu.iupui.proteomics.base.theosp.BpSpec;
import edu.iupui.proteomics.spec.peak.PeakTolerance;
import edu.iupui.proteomics.spec.peak.PositionComparator;
import edu.iupui.proteomics.spec.sp.Ms;

public class DenovoAlign {
	public static Logger logger = Logger.getLogger(DenovoAlign.class);

	private DenovoMng mng;

	/* top-down spectra */
	private Ms<TdDenovoPeak> tdSps[];
	private double precMassMinusWater;

	/* bottom-up de novo resulting sequences */
	private BottomUpSeq allBuSeqs[];
	private BottomUpSeq allEnzymeBuSeqs[][];
	private BottomUpSeq confBuSeqs[];
	private BottomUpSeq enzymeBuSeqs[][];

	/* reference sequence */
	private BpSpec refSeqSpec;

	public static int peptideCount = 0;

	private ArrayList<BottomUpSeq> seleSeqList;

	public DenovoAlign(DenovoMng mng, Ms<TdDenovoPeak> tdPrmSps[],
			BottomUpSeq buSeqs[], BpSpec refSeqSpec) throws Exception {
		this.mng = mng;
		this.tdSps = tdPrmSps;

		this.refSeqSpec = refSeqSpec;
		double refMasses[] = refSeqSpec.getExtBMasses();
		for (int i = 0; i < refMasses.length - 2; i++) {
			System.out.println(refSeqSpec.getResSeq().getAcidString().charAt(i)
					+ " " + refMasses[i + 1]);
		}
		/* initialize enzymeSeqs */
		this.allBuSeqs = buSeqs;
		ArrayList<BottomUpSeq> confSeqList = new ArrayList<BottomUpSeq>();
		for (int i = 0; i < allBuSeqs.length; i++) {
			if (allBuSeqs[i].getDenovoScore() >= mng.confBuSeqThresh) {
				// allBuSeqs[i].output();
				// System.out.println("de novo score " +
				// allBuSeqs[i].getDenovoScore());
				confSeqList.add(allBuSeqs[i]);
			}
		}
		confBuSeqs = confSeqList.toArray(new BottomUpSeq[0]);
		EnumEnzyme enzymes[] = EnumEnzyme.getEnzymeList();
		allEnzymeBuSeqs = new BottomUpSeq[enzymes.length][];
		enzymeBuSeqs = new BottomUpSeq[enzymes.length][];
		for (int i = 0; i < enzymes.length; i++) {
			ArrayList<BottomUpSeq> seqList = new ArrayList<BottomUpSeq>();
			for (int j = 0; j < confBuSeqs.length; j++) {
				if (confBuSeqs[j].getBpSpec().getResSeq().getName()
						.equals(enzymes[i].getName())) {
					seqList.add(confBuSeqs[j]);
				}
			}
			enzymeBuSeqs[i] = seqList.toArray(new BottomUpSeq[0]);

			ArrayList<BottomUpSeq> allSeqList = new ArrayList<BottomUpSeq>();
			for (int j = 0; j < allBuSeqs.length; j++) {
				if (allBuSeqs[j].getBpSpec().getResSeq().getName()
						.equals(enzymes[i].getName())) {
					allSeqList.add(allBuSeqs[j]);
				}
			}
			allEnzymeBuSeqs[i] = allSeqList.toArray(new BottomUpSeq[0]);
		}

		/* get average precursor mass */
		double sumMass = 0;
		for (int i = 0; i < tdSps.length; i++) {
			sumMass += tdSps[i].getHeader().getPrecMonoMassMinusWater();
		}
		precMassMinusWater = sumMass / tdSps.length;
	}

	private ArrayList<AlignPeak> initAlignPeaks() {
		// for (int i = 0; i < tdSps.length; i++) {
		// System.out.println("top down sp " + i);
		// ArrayList<AlignPeak> spList = new ArrayList<AlignPeak>();
		// for (int j = 0; j < tdSps[i].size(); j++) {
		// spList.add(new AlignPeak(tdSps[i].get(j)));
		// }
		// compareRef(spList, spList.size());
		// }
		logger.info("Combine top-down spectra");
		ArrayList<TdDenovoPeak> tdPeakList = new ArrayList<TdDenovoPeak>();
		for (int i = 0; i < tdSps.length; i++) {
			for (int j = 0; j < tdSps[i].size(); j++) {
				tdPeakList.add(tdSps[i].get(j));
			}
		}
		logger.info("combined list length " + tdPeakList.size());
		ArrayList<AlignPeak> peakList = new ArrayList<AlignPeak>();
		for (int i = 0; i < tdPeakList.size(); i++) {
			peakList.add(new AlignPeak(tdPeakList.get(i)));
		}
		return peakList;
	}

	public ArrayList<AlignPeak> removeSimilarPeak(
			ArrayList<AlignPeak> peakList, Comparator<AlignPeak> comp) {
		Collections.sort(peakList, comp);
		ArrayList<AlignPeak> newList = new ArrayList<AlignPeak>();
		for (int i = 0; i < peakList.size(); i++) {
			boolean found = false;
			for (int j = 0; j < newList.size(); j++) {
				double diff = Math.abs(newList.get(j).getPosition()
						- peakList.get(i).getPosition());
				if (diff <= mng.peakMergeTolerance) {
					found = true;
					newList.get(j).mergePeak(peakList.get(i));
					break;
				}
			}
			if (!found) {
				newList.add(peakList.get(i));
				logger.trace("new peak " + newList.size() + " "
						+ peakList.get(i).getPosition() + " "
						+ peakList.get(i).getIntensity());
			}
		}
		Collections.sort(newList, comp);
		return newList;
	}

	public void analyzePeakList(ArrayList<AlignPeak> peakList) {
		for (int support = 1; support <= 10; support++) {
			ArrayList<AlignPeak> newList = new ArrayList<AlignPeak>();
			for (int i = 0; i < peakList.size(); i++) {
				if (peakList.get(i).getTdPeakLen() >= support) {
					newList.add(peakList.get(i));
				}
			}
			System.out.println("Analysis of top down spectra, support number "
					+ support);
			compareRef(newList, newList.size());
		}
	}

	public ArrayList<AlignPeak> removeLowSupportPeak(
			ArrayList<AlignPeak> peakList) {
		ArrayList<AlignPeak> newList = new ArrayList<AlignPeak>();
		for (int i = 0; i < peakList.size(); i++) {
			if (peakList.get(i).getTdPeakLen() >= mng.minSupport) {
				newList.add(peakList.get(i));
			}
		}
		return newList;
	}

	public ArrayList<AlignPeak> removeTdOnlyPeak(ArrayList<AlignPeak> peakList) {
		ArrayList<AlignPeak> newList = new ArrayList<AlignPeak>();
		for (int i = 0; i < peakList.size(); i++) {
			if (peakList.get(i).getBuSeqLen() > 0) {
				newList.add(peakList.get(i));
			}
		}
		return newList;
	}

	public ArrayList<AlignPeak> removeNeighbor(ArrayList<AlignPeak> peakList,
			Comparator<AlignPeak> comp) {
		Collections.sort(peakList, comp);
		boolean keep[] = new boolean[peakList.size()];
		Arrays.fill(keep, true);
		for (int i = 0; i < peakList.size(); i++) {
			if (keep[i]) {
				for (int j = i + 1; j < peakList.size(); j++) {
					double diff = Math.abs(peakList.get(i).getPosition()
							- peakList.get(j).getPosition());
					double error = Math.abs(diff
							- MassConstant.getIsotopeMass());
					if (error < mng.peakTolerance.compRelaxErrorTole(peakList
							.get(i).getPosition(), peakList.get(j)
							.getPosition())) {
						keep[j] = false;
					}
				}
			}
		}
		ArrayList<AlignPeak> newList = new ArrayList<AlignPeak>();
		for (int i = 0; i < peakList.size(); i++) {
			if (keep[i]) {
				newList.add(peakList.get(i));
			}
		}
		logger.info("After removing neighbor size is " + newList.size());
		return newList;
	}

	public ArrayList<AlignPeak> removePeakByAlignScore(
			ArrayList<AlignPeak> peakList, int thresh) {
		ArrayList<AlignPeak> newList = new ArrayList<AlignPeak>();
		boolean remove[] = new boolean[peakList.size()];
		for (int i = 0; i < peakList.size(); i++) {
			if (peakList.get(i).getAlignScore() <= thresh) {
				// && peakList.get(i+1).getAlignScore() <= mng.minAlignScore) {
				remove[i] = true;
				// remove[i+1] = true;
			}
		}
		for (int i = 0; i < peakList.size(); i++) {
			if (!remove[i]) {
				newList.add(peakList.get(i));
			}
		}
		return newList;
	}

	private double[] convertPeakList(ArrayList<AlignPeak> peakList, int num) {
		double masses[] = new double[num];
		for (int i = 0; i < num; i++) {
			masses[i] = peakList.get(i).getPosition();
		}
		return masses;
	}

	private ArrayList<AlignPeak> alignBuSeqs(ArrayList<AlignPeak> peakList) {
		ArrayList<AlignPeak> newList = new ArrayList<AlignPeak>();
		newList.addAll(peakList);
		int count = 0;
		System.out.println("Number of bottom-up " + confBuSeqs.length);
		for (int i = 0; i < confBuSeqs.length; i++) {
			confBuSeqs[i].findBestShift(
					convertPeakList(peakList, peakList.size()),
					mng.peakTolerance);
			if (confBuSeqs[i].getBestScore() >= 5
					&& confBuSeqs[i].getBestShift() >= 20000) {
				confBuSeqs[i].output();
				compareRef(confBuSeqs[i].getShiftedPrms(),
						confBuSeqs[i].getShiftedPrms().length);
			}
			// System.out.println("score " + confBuSeqs[i].getBestScore());
			if (confBuSeqs[i].getBestScore() >= mng.minBuScore) {

				count++;
				double masses[] = confBuSeqs[i].getShiftedPrms();
				for (int j = 0; j < masses.length; j++) {
					if (masses[j] > mng.minMass
							&& masses[j] < precMassMinusWater - mng.minMass) {
						newList.add(new AlignPeak(masses[j], confBuSeqs[i]));
					}
				}
			}
		}
		System.out
				.println("Number of bottom-up aligned to top-down (threshold) "
						+ mng.minBuScore + " is " + count);
		return newList;
	}

	private ArrayList<AlignPeak> removeBuPeakList(ArrayList<AlignPeak> peakList) {
		ArrayList<AlignPeak> newList = new ArrayList<AlignPeak>();
		for (int i = 0; i < peakList.size(); i++) {
			AlignPeak peak = peakList.get(i);
			if (peak.getTdPeakLen() != 0) {
				peak.removeAllBuSeqs();
				newList.add(peak);
			}
		}
		return newList;
	}

	private ArrayList<AlignPeak> coverAlignBuSeqs(ArrayList<AlignPeak> peakList) {
		@SuppressWarnings("unchecked")
		ArrayList<BottomUpSeq> seqList[] = new ArrayList[(int) Math
				.round(precMassMinusWater / 100) + 1];
		for (int i = 0; i < seqList.length; i++) {
			seqList[i] = new ArrayList<BottomUpSeq>();
		}

		BottomUpSeq uniqBuSeqs[] = BottomUpSeq.removeSimilarMass(confBuSeqs);
		logger.debug("Add buseqs. buseqs number " + uniqBuSeqs.length);
		for (int i = 0; i < uniqBuSeqs.length; i++) {
			uniqBuSeqs[i].findEndBestShift(
					convertPeakList(peakList, peakList.size()),
					mng.peakTolerance);
			uniqBuSeqs[i].output();
			double masses[] = uniqBuSeqs[i].getShiftedPrms();
			if (masses.length > 0) {
				double left = masses[0];
				double right = masses[masses.length - 1];
				int lpos = (int) Math.round(left / 100);
				int rpos = (int) Math.round(right / 100);
				for (int j = lpos; j <= rpos; j++) {
					if (j > 0 && j < seqList.length) {
						seqList[j].add(uniqBuSeqs[i]);
					}
				}
			}
		}

		seleSeqList = new ArrayList<BottomUpSeq>();
		for (int i = 0; i < seqList.length; i++) {
			Collections.sort(seqList[i], new NormScoreComparator());
			for (int j = 0; j < mng.seleSeqNum; j++) {
				if (j < seqList[i].size()
						&& seqList[i].get(j).getNormalizedScore() >= 0.4
						&& !seleSeqList.contains(seqList[i].get(j))) {
					seleSeqList.add(seqList[i].get(j));
					System.out.println("Add new ");
					seqList[i].get(j).output();
					compareRef(seqList[i].get(j).getShiftedTagPrms(),
							seqList[i].get(j).getShiftedTagPrms().length);
				}
			}
		}

		peakList = removeBuPeakList(peakList);
		for (int i = 0; i < seleSeqList.size(); i++) {
			double masses[] = seleSeqList.get(i).getShiftedTagPrms();
			for (int j = 0; j < masses.length; j++) {
				if (masses[j] > mng.minMass
						&& masses[j] < precMassMinusWater - mng.minMass) {
					peakList.add(new AlignPeak(masses[j], seleSeqList.get(i)));
				}
			}
		}
		return peakList;
	}

	public ArrayList<AlignPeak> keepTopPeaks(ArrayList<AlignPeak> peakList) {
		ArrayList<AlignPeak> tdList = new ArrayList<AlignPeak>();
		ArrayList<AlignPeak> buList = new ArrayList<AlignPeak>();
		for (int i = 0; i < peakList.size(); i++) {
			if (peakList.get(i).getTdPeakLen() > 0) {
				tdList.add(peakList.get(i));
			} else {
				buList.add(peakList.get(i));
			}
		}
		ArrayList<AlignPeak> newList = new ArrayList<AlignPeak>();
		int count = 0;
		for (int i = 0; i < tdList.size(); i++) {
			if (tdList.get(i).getBuSeqLen() >= 1) {
				newList.add(tdList.get(i));
				count++;
			}
		}
		logger.info(count + " td prms are kept.");
		count = 0;
		for (int i = 0; i < buList.size(); i++) {
			if (buList.get(i).getPosition() <= precMassMinusWater - mng.minMass
					&& buList.get(i).getPosition() >= mng.minMass
					&& buList.get(i).getBuSeqLen() >= 2) {
				newList.add(buList.get(i));
				count++;
			}
		}
		logger.info(count + " bu prms are kept.");
		return newList;
	}

	public ArrayList<AlignPeak> removeTdComplementary(
			ArrayList<AlignPeak> peakList) {
		ArrayList<AlignPeak> tdList = new ArrayList<AlignPeak>();
		ArrayList<AlignPeak> buList = new ArrayList<AlignPeak>();
		for (int i = 0; i < peakList.size(); i++) {
			if (peakList.get(i).getTdPeakLen() > 0) {
				tdList.add(peakList.get(i));
			} else {
				buList.add(peakList.get(i));
			}
		}
		boolean keep[] = new boolean[tdList.size()];
		Arrays.fill(keep, true);
		for (int i = 0; i < tdList.size(); i++) {
			if (keep[i]) {
				for (int j = i + 1; j < tdList.size(); j++) {
					double sum = tdList.get(i).getPosition()
							+ tdList.get(j).getPosition();
					double error = Math.abs(sum - precMassMinusWater
							- MassConstant.getWaterMass());
					if (error < mng.peakTolerance
							.compStrictErrorTole(precMassMinusWater)) {
						keep[j] = false;
					}
				}
			}
		}
		ArrayList<AlignPeak> newList = new ArrayList<AlignPeak>();
		for (int i = 0; i < tdList.size(); i++) {
			if (keep[i]) {
				newList.add(tdList.get(i));
			}
		}
		for (int i = 0; i < buList.size(); i++) {
			newList.add(buList.get(i));
		}
		logger.info("After removing complementary size is " + newList.size());
		return newList;
	}

	private boolean findMass(ArrayList<AlignPeak> peakList, double mass,
			PeakTolerance peakTolerance) {
		double tolerance = peakTolerance.compStrictErrorTole(mass);
		for (int i = 0; i < peakList.size(); i++) {
			double diff = Math.abs(mass - peakList.get(i).getPosition());
			if (diff <= tolerance) {
				return true;
			}
		}
		return false;
	}

	private AlignPeak findMatchPeak(double mass, ArrayList<AlignPeak> peakList,
			PeakTolerance peakTolerance) throws Exception {
		double tolerance = peakTolerance.compStrictErrorTole(mass);
		for (int i = 0; i < peakList.size(); i++) {
			double diff = Math.abs(mass - peakList.get(i).getPosition());
			if (diff <= tolerance) {
				return peakList.get(i);
			}
		}
		return null;
	}

	private int compSuppNum(double mass, ArrayList<AlignPeak> peakList,
			PeakTolerance peakTolerance) throws Exception {
		ResList resList = ResListFactory.getSystemInstance(mng.buResFileName);
		boolean leftSupport = false;
		boolean rightSupport = false;
		for (int i = 0; i < 20; i++) {
			double resMass = resList.get(i).getMass();
			double leftMass = mass - resMass;
			if (findMass(peakList, leftMass, peakTolerance)) {
				leftSupport = true;
			}
			double rightMass = mass + resMass;
			if (findMass(peakList, rightMass, peakTolerance)) {
				rightSupport = true;
			}
		}
		if (leftSupport && rightSupport) {
			return 2;
		} else if (leftSupport || rightSupport) {
			return 1;
		} else {
			return 0;
		}
	}

	public ArrayList<AlignPeak> adjustTdPeak(ArrayList<AlignPeak> peakList,
			ArrayList<AlignPeak> tdPeakList) throws Exception {
		for (int i = 0; i < tdPeakList.size(); i++) {
			TdDenovoPeak tdPeak = tdPeakList.get(i).getTdPeakList().get(0);
			double curMass = tdPeak.getPosition();
			AlignPeak matchPeak = findMatchPeak(curMass, peakList,
					mng.peakTolerance);
			if (matchPeak != null) {
				continue;
			}
			double leftMass = curMass - MassConstant.getIsotopeMass();
			double rightMass = curMass + MassConstant.getIsotopeMass();
			AlignPeak leftPeak = findMatchPeak(leftMass, peakList,
					mng.peakTolerance);
			AlignPeak rightPeak = findMatchPeak(rightMass, peakList,
					mng.peakTolerance);
			int curSuppNum = compSuppNum(curMass, peakList, mng.peakTolerance);
			if (leftPeak == null && rightPeak == null) {
				continue;
			}
			if (leftPeak != null && rightPeak != null) {
				continue;
			}
			AlignPeak peak;
			double newMass;
			if (leftPeak != null) {
				peak = leftPeak;
				newMass = leftMass;
			} else {
				peak = rightPeak;
				newMass = rightMass;
			}

			logger.debug("Mass " + curMass + " " + leftPeak + " " + curSuppNum
					+ " " + rightPeak + " new mass " + newMass);
			ArrayList<TdDenovoPeak> list = tdPeakList.get(i).getTdPeakList();
			for (int j = 0; j < list.size(); j++) {
				list.get(j).setPosition(newMass);
			}
			peak.mergePeak(tdPeakList.get(i));
		}
		return peakList;
	}

	public ArrayList<AlignPeak> adjustTdPeak(ArrayList<AlignPeak> peakList)
			throws Exception {
		ArrayList<AlignPeak> tdList = new ArrayList<AlignPeak>();
		for (int i = 0; i < peakList.size(); i++) {
			if (peakList.get(i).getTdPeakLen() > 0) {
				tdList.add(peakList.get(i));
			}
		}
		for (int i = 0; i < tdList.size(); i++) {
			TdDenovoPeak tdPeak = tdList.get(i).getTdPeakList().get(0);
			double curMass = tdPeak.getPosition();
			double leftMass = curMass - MassConstant.getIsotopeMass();
			double rightMass = curMass + MassConstant.getIsotopeMass();
			int curSuppNum = compSuppNum(curMass, peakList, mng.peakTolerance);
			int leftSuppNum = compSuppNum(leftMass, peakList, mng.peakTolerance);
			int rightSuppNum = compSuppNum(rightMass, peakList,
					mng.peakTolerance);

			double newMass = curMass;
			if (leftSuppNum - curSuppNum == 2) {
				newMass = leftMass;
			} else if (rightSuppNum - curSuppNum == 2) {
				newMass = rightMass;
			}
			logger.debug("Mass " + curMass + " " + leftSuppNum + " "
					+ curSuppNum + " " + rightSuppNum + " new mass " + newMass);
			tdList.get(i).setPosition(newMass);
			ArrayList<TdDenovoPeak> list = tdList.get(i).getTdPeakList();
			for (int j = 0; j < list.size(); j++) {
				list.get(j).setPosition(newMass);
			}
		}
		return peakList;
	}

	private Interval getGap(ArrayList<AlignPeak> peakList, double start,
			double gapLen) {
		int leftIdx = -1;
		int rightIdx = -1;
		double left = start;
		double right = start;
		for (int i = 0; i < peakList.size(); i++) {
			double pos = peakList.get(i).getPosition();
			if (pos > left) {
				if (pos - left <= gapLen) {
					left = pos;
					leftIdx = i;
				} else {
					right = pos;
					rightIdx = i;
					break;
				}
			}
		}
		if (right - left <= gapLen) {
			return null;
		} else {
			if (leftIdx > 0) {
				left = peakList.get(leftIdx).getPosition();
			}
			if (rightIdx < peakList.size() - 1) {
				right = peakList.get(rightIdx).getPosition();
			}
		}
		return new Interval(left, right);
	}

	// private Interval getLongGap(ArrayList<AlignPeak> peakList, double start,
	// double gapLen) {
	// if (start < 6994) {
	// return new Interval (6894.556163125, 10610.340808460938);
	// }
	// else {
	// return null;
	// }
	// }

	public ArrayList<AlignPeak> getAddPeakList(
			ArrayList<AlignPeak> oriPeakList, ArrayList<AlignPeak> peakList,
			double min, double max) {
		ArrayList<AlignPeak> addList = new ArrayList<AlignPeak>();
		int count[] = new int[(int) Math.round(max / 100) + 1];
		for (int i = 0; i < oriPeakList.size(); i++) {
			double pos = oriPeakList.get(i).getPosition();
			if (pos < min || pos > max) {
				continue;
			}
			boolean found = false;
			double error = mng.peakTolerance
					.compStrictErrorTole(precMassMinusWater);
			for (int j = 0; j < peakList.size(); j++) {
				double bp = peakList.get(j).getPosition();
				double diff = Math.abs(pos - bp);
				if (diff <= error) {
					found = true;
					break;
				}
				if (Math.abs(diff - MassConstant.getIsotopeMass()) <= error) {
					found = true;
					break;
				}
				double comp_bp = precMassMinusWater
						+ MassConstant.getWaterMass() - bp;
				diff = Math.abs(pos - comp_bp);
				if (diff <= error) {
					found = true;
					break;
				}
				if (Math.abs(diff - MassConstant.getIsotopeMass()) <= error) {
					found = true;
					break;
				}

			}
			if (!found) {
				int index = (int) Math
						.round(oriPeakList.get(i).getPosition() / 100);
				if (count[index] < 3) {
					// oriPeakList.get(i).setTdConf(AlignPeak.LOW_TD_CONF);
					addList.add(oriPeakList.get(i));
					count[index]++;
				}
			}
		}
		return addList;
	}

	public ArrayList<AlignPeak> getNeighborList(ArrayList<AlignPeak> peakList,
			double min, double max, double expand) {
		ArrayList<AlignPeak> neighborList = new ArrayList<AlignPeak>();
		for (int i = 0; i < peakList.size(); i++) {
			double pos = peakList.get(i).getPosition();
			if (pos <= min && pos > min - expand) {
				AlignPeak peak = peakList.get(i);
				peak.setLeftAnchor(true);
				peak.setRightAnchor(false);
				neighborList.add(peak);
			}

			if (pos >= max && pos < max + expand) {
				AlignPeak peak = peakList.get(i);
				peak.setLeftAnchor(false);
				peak.setRightAnchor(true);
				neighborList.add(peak);
			}
		}
		return neighborList;
	}

	private ArrayList<AlignPeak> fillGap(ArrayList<AlignPeak> oriPeakList,
			ArrayList<AlignPeak> peakList, Interval gap, double flankLen)
			throws Exception {
		ArrayList<AlignPeak> tdPeakList = getAddPeakList(oriPeakList, peakList,
				gap.getLeft(), gap.getRight());
		Collections.sort(tdPeakList, new PositionComparator());
		logger.debug("Top-down peak list " + tdPeakList.size());
		ArrayList<AlignPeak> neighborList = getNeighborList(peakList,
				gap.getLeft(), gap.getRight(), flankLen);
		logger.debug("Neighbor peak List " + neighborList.size());
		ArrayList<AlignPeak> curPeakList = new ArrayList<AlignPeak>();
		curPeakList.addAll(tdPeakList);
		curPeakList.addAll(neighborList);
		Interval bound = new Interval(gap.getLeft() - flankLen, gap.getRight()
				+ flankLen);
		boolean stop = false;
		int count = 0;
		peptideCount = 0;
		while (!stop) {
			System.out.println("EM round " + (count + 1) + " bound left "
					+ bound.getLeft() + " right " + bound.getRight());
			stop = EM(curPeakList, gap, bound);
			count++;
			if (count == 60) {
				stop = true;
			}
		}
		System.out.println("Peptide count " + peptideCount);
		ArrayList<AlignPeak> newList = new ArrayList<AlignPeak>();
		newList.addAll(peakList);
		for (int i = 0; i < curPeakList.size(); i++) {
			AlignPeak peak = curPeakList.get(i);
			if (peak.getPosition() > gap.getLeft()
					&& peak.getPosition() < gap.getRight()) {
				newList.add(peak);
			}
		}
		Collections.sort(newList, new PositionComparator());
		return newList;
	}

	private ArrayList<AlignPeak> fillLongGap(ArrayList<AlignPeak> oriPeakList,
			ArrayList<AlignPeak> peakList, double leftBound, double rightBound,
			double flankLen) throws Exception {
		BottomUpSeq bestSeq = null;
		double bestScore = 0;
		EnumEnzyme enzymes[] = EnumEnzyme.getEnzymeList();
		for (int i = 0; i < enzymes.length; i++) {
			EnumEnzyme enzyme = enzymes[i];
			logger.debug("Enzyme " + enzyme.getName());
			int id = enzyme.getId();
			ArrayList<Double> leftMasses = getLongGapMasses(peakList, enzyme,
					leftBound - flankLen, leftBound);
			ArrayList<Double> rightMasses = getLongGapMasses(peakList, enzyme,
					rightBound, rightBound + flankLen);
			rightMasses.add(rightBound);

			for (int j = 0; j < leftMasses.size(); j++) {
				for (int k = 0; k < rightMasses.size(); k++) {
					double theoMass = rightMasses.get(k) - leftMasses.get(j);
					logger.trace("theoretical mass " + theoMass + " left mass "
							+ leftMasses.get(j) + " right mass "
							+ rightMasses.get(k));
					for (int l = 0; l < allEnzymeBuSeqs[id].length; l++) {
						BottomUpSeq seq = allEnzymeBuSeqs[id][l];
						double pepMass = seq.getSimilarResSumMass(theoMass);
						double error = Math.abs(pepMass - theoMass);
						if (error <= mng.longGapPrecTolerance) {
							logger.trace("Found bottom-up sequence "
									+ seq.getBpSpec().getResSeq()
											.getAcidString() + " "
									+ seq.getDenovoScore());
							if (seq.getDenovoScore() > bestScore) {
								bestScore = seq.getDenovoScore();
								bestSeq = seq;
								bestSeq.setBestShift(leftMasses.get(j));
							}
						}

					}
				}
			}
		}

		ArrayList<AlignPeak> tdPeakList = getAddPeakList(oriPeakList, peakList,
				leftBound, rightBound);
		for (int i = 0; i < tdPeakList.size(); i++) {
			logger.debug("Add ori peak " + tdPeakList.get(i).getPosition());
		}
		peakList.addAll(tdPeakList);

		if (bestSeq != null) {
			double prms[] = bestSeq.getCysShiftedPrms(22);
			logger.debug("Best sequence "
					+ bestSeq.getBpSpec().getResSeq().getAcidString()
					+ " shift " + bestSeq.getBestShift());
			for (int j = 0; j < prms.length; j++) {
				if (prms[j] > mng.minMass
						&& prms[j] < precMassMinusWater - mng.minMass) {
					logger.debug("Add peak " + prms[j]);
					AlignPeak peak = new AlignPeak(prms[j], bestSeq);
					peakList.add(peak);
				}
			}
		}
		return peakList;
	}

	public void process() throws Exception {
		/* step 1 Get peak list from top-down spectra */
		ArrayList<AlignPeak> tdPeakList = initAlignPeaks();
		tdPeakList = removeSimilarPeak(tdPeakList, new TdPeakComparator());
		logger.info("list length after removing similar masses "
				+ tdPeakList.size());
		tdPeakList = removeNeighbor(tdPeakList, new TdPeakComparator());
		analyzePeakList(tdPeakList);
		tdPeakList = removeLowSupportPeak(tdPeakList);
		logger.info("list length after removing low support peaks "
				+ tdPeakList.size());
		compareRef(tdPeakList, tdPeakList.size());
		Collections.sort(tdPeakList, new PositionComparator());
		outputPeakList(tdPeakList);

		/* step 2 align bottom-up to peak list */
		logger.info("Add buseqs. buseqs number " + confBuSeqs.length);
		ArrayList<AlignPeak> peakList = alignBuSeqs(tdPeakList);
		compareRef(peakList, peakList.size());
		peakList = removeSimilarPeak(peakList, new TdBuPeakComparator());
		peakList = keepTopPeaks(peakList);
		peakList = removeTdComplementary(peakList);
		peakList = removeNeighbor(peakList, new BuTdComparator());
		compareRef(peakList, peakList.size());
		Collections.sort(peakList, new PositionComparator());
		outputPeakList(peakList);

		// step 3 fill gap phase 1 
		// import sort td peaks based on tdPeakComparator 
		Collections.sort(tdPeakList, new TdPeakComparator());
		logger.info("Start filling gaps.");
		double start = 0;
		Interval gap;
		while ((gap = getGap(peakList, start, mng.gapLen)) != null) {
			logger.debug("Gap found " + gap.getLeft() + " to " + gap.getRight());
			peakList = fillGap(tdPeakList, peakList, gap, mng.flankLen);
			compareRef(peakList, peakList.size());
			start = gap.getRight();
		}
		compareRef(peakList, peakList.size());
		// peakList = adjustTdPeak(peakList);
		compareRef(peakList, peakList.size());
		Collections.sort(peakList, new PositionComparator());
		outputPeakList(peakList);
		// repeat step 2 
		peakList = alignBuSeqs(peakList);
		compareRef(peakList, peakList.size());
		peakList = removeSimilarPeak(peakList, new TdBuPeakComparator());
		peakList = keepTopPeaks(peakList);
		peakList = removeTdComplementary(peakList);
		peakList = removeNeighbor(peakList, new BuTdComparator());
		compareRef(peakList, peakList.size());
		Collections.sort(peakList, new PositionComparator());
		outputPeakList(peakList);

		// step 4 fill gap phase 2
		start = 0;
		System.out.println("ORI TD PEAK LIST");
		outputPeakList(tdPeakList);
		while ((gap = getGap(peakList, start, mng.longGapLen)) != null) {
			logger.debug("Find long gap " + gap.getLeft() + " "
					+ gap.getRight());
			fillLongGap(tdPeakList, peakList, gap.getLeft(), gap.getRight(),
					mng.flankLen);
			start = gap.getRight();
		}
		Collections.sort(peakList, new PositionComparator());
		compareRef(peakList, peakList.size());
		// outputPeakList(peakList);

		peakList = coverAlignBuSeqs(peakList);
		compareRef(peakList, peakList.size());
		Collections.sort(peakList, new PositionComparator());
		peakList = removeSimilarPeak(peakList, new TdBuPeakComparator());
		peakList = removeTdComplementary(peakList);
		peakList = removeNeighbor(peakList, new BuTdComparator());
		compareRef(peakList, peakList.size());
		Collections.sort(peakList, new PositionComparator());
		outputPeakList(peakList);

		peakList = adjustTdPeak(peakList);
		peakList = adjustTdPeak(peakList, tdPeakList);
		peakList = removeTdOnlyPeak(peakList);
		System.out.println("Adjust td peaks");
		compareRef(peakList, peakList.size());
		outputPeakList(peakList);
		

		/* Start find path */
		FindPath findPath = new FindPath(peakList, precMassMinusWater, mng);
		Res residues[] = findPath.process();
		correctKandQ(residues);
		ArrayList<AlignPeak> finalPeakList = findPath.getResultPeakList();
		compareRef(finalPeakList, finalPeakList.size());
		outputPeakList(finalPeakList);

		System.out.println("Path result");
		for (int i = 0; i < residues.length; i++) {
			System.out.print(residues[i].getAcid().getOneLetter());
		}
		System.out.println();
		peakList = findPath.getResultPeakList();
		Collections.sort(peakList, new PositionComparator());
		compareRef(peakList, peakList.size());
		outputPeakList(peakList);
		peakList = removePeakByAlignScore(peakList, 1);
		compareRef(peakList, peakList.size());
		outputPeakList(peakList);
	}

	private void correctKandQ(Res residues[]) throws Exception {
		double mass = 0;
		ResList resList = ResListFactory.getSystemInstance(mng.buResFileName);
		for (int i = 0; i < residues.length; i++) {

			String letter = residues[i].getAcid().getOneLetter();
			if (letter.equals("K") || letter.equals("Q")) {
				int k_count = 0;
				int q_count = 0;
				for (int j = 0; j < seleSeqList.size(); j++) {
					double buMasses[] = seleSeqList.get(j).getShiftedPrms();
					for (int k = 0; k < buMasses.length - 1; k++) {
						if (Math.abs(buMasses[k] - mass) <= 0.2) {
							int pos = k;
							if (seleSeqList.get(j).getNGap() == 0.0) {
								pos = pos - 1;
							}
							if (pos >= 0
									&& pos < seleSeqList.get(j).getBpSpec()
											.getResSeq().getLen()) {
								String buLetter = seleSeqList.get(j)
										.getBpSpec().getResSeq().getRes(k)
										.getAcid().getOneLetter();
								if (buLetter.equals("K")) {
									k_count++;
								} else if (buLetter.equals("Q")) {
									q_count++;
								}
							}
						}
					}
				}
				System.out.println("pos " + i + " k count " + k_count
						+ " q count " + q_count);
				if (k_count > q_count) {
					residues[i] = resList.getResByAcid("K");
				} else if (q_count > k_count) {
					residues[i] = resList.getResByAcid("Q");
				}
			}
			mass = mass + residues[i].getMass();
		}
	}

	/********************** support functions ********************/
	private void compareRef(ArrayList<AlignPeak> peakList, int peakNum) {
		ArrayList<AlignPeak> seleList = new ArrayList<AlignPeak>();
		for (int i = 0; i < peakNum; i++) {
			seleList.add(peakList.get(i));
			peakList.get(i).setLeftMatch(false);
			peakList.get(i).setExactMatch(false);
			peakList.get(i).setRightMatch(false);
		}
		Collections.sort(seleList, new PositionComparator());
		String seq = refSeqSpec.getResSeq().getAcidString();
		double refPrms[] = refSeqSpec.getExtBMasses();
		boolean isLeftMatch[] = new boolean[refPrms.length];
		boolean isExactMatch[] = new boolean[refPrms.length];
		boolean isRightMatch[] = new boolean[refPrms.length];
		int count = 0;
		for (int i = 0; i < refPrms.length; i++) {
			for (int j = 0; j < peakNum; j++) {
				double tolerance = mng.peakTolerance.compRelaxErrorTole(
						refPrms[i], precMassMinusWater);
				if (Math.abs(peakList.get(j).getPosition() - refPrms[i]) <= tolerance) {
					peakList.get(j).setExactMatch(true);
					isExactMatch[i] = true;
					count++;
				}
				if (Math.abs(peakList.get(j).getPosition() - refPrms[i]
						- 1.00235) <= tolerance) {
					peakList.get(j).setRightMatch(true);
					isRightMatch[i] = true;
					count++;
				}
				if (Math.abs(peakList.get(j).getPosition() - refPrms[i]
						+ 1.00235) <= tolerance) {
					peakList.get(j).setLeftMatch(true);
					isLeftMatch[i] = true;
					count++;
				}
			}
		}

		System.out.println("Total peak " + peakNum + " matched peak: " + count
				+ " percentage " + (count / (float) peakNum));
		String anno = "";
		int num = 0;
		for (int i = 0; i < seq.length(); i++) {
			anno = anno + seq.charAt(i);
			if (isExactMatch[i + 1]) {
				anno = anno + "|";
				num++;
			} else {
				anno = anno + " ";
			}
		}
		System.out.println("Exact match  " + num + " anno " + anno);
		anno = "";
		num = 0;
		for (int i = 0; i < seq.length(); i++) {
			anno = anno + seq.charAt(i);
			if (isLeftMatch[i + 1]) {
				anno = anno + "|";
				num++;
			} else {
				anno = anno + " ";
			}
		}
		System.out.println("Left match  " + num + " anno " + anno);
		anno = "";
		num = 0;
		for (int i = 0; i < seq.length(); i++) {
			anno = anno + seq.charAt(i);
			if (isRightMatch[i + 1]) {
				anno = anno + "|";
				num++;
			} else {
				anno = anno + " ";
			}
		}
		System.out.println("Right match  " + num + " anno " + anno);
		anno = "";
		num = 0;
		for (int i = 0; i < seq.length(); i++) {
			anno = anno + seq.charAt(i);
			if (isExactMatch[i + 1] || isLeftMatch[i + 1]
					|| isRightMatch[i + 1]) {
				anno = anno + "|";
				num++;
			} else {
				anno = anno + " ";
			}
		}
		System.out.println("Combine match " + num + " anno " + anno);
	}

	public void outputBuSeq(BottomUpSeq seq, ArrayList<AlignPeak> peakList) {
		System.out.println("shift " + seq.getBestShift());
		double masses[] = seq.getShiftedPrms();
		for (int i = 0; i < masses.length - 1; i++) {
			System.out.println(i + " "
					+ seq.getBpSpec().getResSeq().getAcidString().charAt(i)
					+ " " + masses[i + 1]);
		}
		double tdmasses[] = convertPeakList(peakList, peakList.size());
		Arrays.sort(tdmasses);
		for (int i = 0; i < tdmasses.length; i++) {
			boolean found = false;
			for (int j = 0; j < masses.length; j++) {
				double diff = Math.abs(tdmasses[i] - masses[j]);
				double error = mng.peakTolerance
						.compStrictErrorTole(tdmasses[i]);
				if (diff <= error) {
					found = true;
					System.out.println("Td mass " + tdmasses[i] + " bu mass "
							+ masses[j]);
					break;
				}
			}
			if (!found) {
				System.out.println("Td mass " + tdmasses[i]);
			}

		}
	}

	public void outputPeakList(ArrayList<AlignPeak> peakList) {
		for (int i = 0; i < peakList.size(); i++) {
			if (!peakList.get(i).isExactMatch()) {
				System.out.print("***");
			}
			System.out.println(i + "\t" + peakList.get(i).getTdPeakLen() + "\t"
					+ peakList.get(i).getBuSeqLen() + "\t"
					+ peakList.get(i).getPosition() + "\t"
					+ (precMassMinusWater - peakList.get(i).getPosition())
					+ "\t" + peakList.get(i).isLeftMatch() + "\t"
					+ peakList.get(i).isExactMatch() + "\t"
					+ peakList.get(i).isRightMatch() + "\t"
					+ peakList.get(i).getAlignScore());
		}
	}

	public boolean findMass(ArrayList<AlignPeak> leftAnchorPeakList,
			double mass, double tolerance) {
		for (int i = 0; i < leftAnchorPeakList.size(); i++) {
			double diff = Math.abs(leftAnchorPeakList.get(i).getPosition()
					- mass);
			if (diff <= tolerance) {
				return true;
			}
		}
		return false;
	}

	public ArrayList<Double> getLongGapMasses(ArrayList<AlignPeak> peakList,
			EnumEnzyme enzyme, double minMass, double maxMass) throws Exception {
		ResList resList = ResListFactory.getSystemInstance(mng.buResFileName);
		logger.debug("Min mass" + minMass + " max mass " + maxMass);
		ArrayList<Double> massList = new ArrayList<Double>();
		for (int i = 1; i < peakList.size(); i++) {
			double mass = peakList.get(i).getPosition();
			if (mass < minMass || mass > maxMass) {
				continue;
			}
			String breakResidue = enzyme.getBreakResidue();
			for (int j = 0; j < breakResidue.length(); j++) {
				String res = breakResidue.substring(j, j + 1);
				double prevMass = mass - resList.getResByAcid(res).getMass();
				if (findMass(peakList, prevMass, mng.fillgGapTolerance)) {
					massList.add(mass);
					continue;
				}
			}
		}
		return massList;
	}

	private Cleavage getCleavage(ArrayList<AlignPeak> peakList, double mass,
			ResList resList, EnumEnzyme enzymes[]) {
		Cleavage cleavage = null;
		for (int enz = 0; enz < enzymes.length; enz++) {
			EnumEnzyme enzyme = enzymes[enz];
			String breakResidue = enzyme.getBreakResidue();
			for (int j = 0; j < breakResidue.length(); j++) {
				String res = breakResidue.substring(j, j + 1);
				double prevMass = mass - resList.getResByAcid(res).getMass();
				if (findMass(peakList, prevMass, mng.fillgGapTolerance)) {
					if (cleavage == null) {
						cleavage = new Cleavage(mass);
					}
					cleavage.addEnzyme(enzyme);
				}
			}
		}
		return cleavage;
	}

	public ArrayList<Cleavage> getLeftCleavages(
			ArrayList<AlignPeak> leftAnchorPeakList, double minMass)
			throws Exception {
		ResList resList = ResListFactory.getSystemInstance(mng.buResFileName);
		ArrayList<Cleavage> cleavageList = new ArrayList<Cleavage>();
		logger.debug("Min mass " + minMass);
		EnumEnzyme enzymes[] = EnumEnzyme.getEnzymeList();

		for (int i = 1; i < leftAnchorPeakList.size(); i++) {
			double mass = leftAnchorPeakList.get(i).getPosition();
			if (mass <= minMass + 1) {
				continue;
			}
			Cleavage cleavage = getCleavage(leftAnchorPeakList, mass, resList,
					enzymes);
			if (cleavage != null) {
				cleavageList.add(cleavage);
			}
		}
		return cleavageList;
	}

	public ArrayList<Cleavage> getRightCleavages(
			ArrayList<AlignPeak> rightAnchorPeakList, double maxMass)
			throws Exception {
		ArrayList<Cleavage> cleavageList = new ArrayList<Cleavage>();

		EnumEnzyme enzymes[] = EnumEnzyme.getEnzymeList();
		ResList resList = ResListFactory.getSystemInstance(mng.buResFileName);
		logger.debug("Max mass " + maxMass);
		for (int i = rightAnchorPeakList.size() - 1; i >= 1; i--) {
			double mass = rightAnchorPeakList.get(i).getPosition();
			if (mass >= maxMass - 1) {
				continue;
			}
			Cleavage cleavage = getCleavage(rightAnchorPeakList, mass, resList,
					enzymes);
			if (cleavage != null) {
				cleavageList.add(cleavage);
				logger.debug("Add cleavage " + cleavage.getMass());
			}
		}
		return cleavageList;
	}

	public ArrayList<AlignPeak> keepLocalTop(ArrayList<AlignPeak> peakList,
			double min, double max) {
		Collections.sort(peakList, new LocalTopComparator());
		ArrayList<AlignPeak> newList = new ArrayList<AlignPeak>();
		int count[] = new int[(int) Math.round(max / 100) + 1];
		int num = 0;
		for (int i = 0; i < peakList.size(); i++) {
			double pos = peakList.get(i).getPosition();
			if (pos < min || pos > max) {
				newList.add(peakList.get(i));
				num++;
			} else {
				int index = (int) Math.round(pos / 100);
				if (count[index] < mng.keepLocalNum) {
					newList.add(peakList.get(i));
					count[index]++;
					num++;
				}
			}
		}
		logger.debug(num + " masses are kept ");
		return newList;
	}

	public boolean EM(ArrayList<AlignPeak> curPeakList, Interval gap,
			Interval bound) throws Exception {
		boolean stop = true;
		ArrayList<AlignPeak> leftAnchorPeakList = new ArrayList<AlignPeak>();
		for (int i = 0; i < curPeakList.size(); i++) {
			if (curPeakList.get(i).isLeftAnchor()) {
				leftAnchorPeakList.add(curPeakList.get(i));
			}
		}
		ArrayList<AlignPeak> rightAnchorPeakList = new ArrayList<AlignPeak>();
		for (int i = 0; i < curPeakList.size(); i++) {
			if (curPeakList.get(i).isRightAnchor()) {
				rightAnchorPeakList.add(curPeakList.get(i));
			}
		}
		Collections.sort(leftAnchorPeakList, new PositionComparator());
		Collections.sort(rightAnchorPeakList, new PositionComparator());
		Collections.sort(curPeakList, new PositionComparator());
		logger.debug("Left anchor size " + leftAnchorPeakList.size()
				+ " right anchor size " + rightAnchorPeakList.size()
				+ " peak size " + curPeakList.size());

		logger.debug("Start EM ");
		stop = true;
		ArrayList<Cleavage> leftCleavage = getLeftCleavages(leftAnchorPeakList,
				bound.getLeft());
		if (leftCleavage.size() > 0) {
			stop = false;
			bound.setLeft(leftCleavage.get(0).getMass());
		}
		ArrayList<AlignPeak> leftNewPeakList = new ArrayList<AlignPeak>();
		double bestScore = -1;
		double bestStart = -1;
		for (int r = 0; r < leftCleavage.size(); r++) {
			logger.debug("Left cleavage mass " + leftCleavage.get(r).getMass());
			double shifts[] = new double[1];
			shifts[0] = leftCleavage.get(r).getMass();
			ArrayList<EnumEnzyme> enzymeList = leftCleavage.get(r)
					.getEnzymeList();
			for (int i = 0; i < enzymeList.size(); i++) {
				EnumEnzyme enzyme = enzymeList.get(i);
				int id = enzyme.getId();
				for (int j = 0; j < enzymeBuSeqs[id].length; j++) {
					enzymeBuSeqs[id][j].findLimitedBestShift(
							convertPeakList(curPeakList, curPeakList.size()),
							shifts, mng.peakTolerance);
					if (enzymeBuSeqs[id][j].getBestScore() > mng.fillGapMinScore) {
						peptideCount++;
						enzymeBuSeqs[id][j].output();
						System.out.println("Best shift "
								+ enzymeBuSeqs[id][j].getBestShift()
								+ " best score "
								+ enzymeBuSeqs[id][j].getBestScore());
					}
					if (enzymeBuSeqs[id][j].getBestScore() > bestScore) {
						bestScore = enzymeBuSeqs[id][j].getBestScore();
						bestStart = leftCleavage.get(r).getMass();
					}
				}
			}
		}
		if (bestScore < mng.fillGapMinScore) {
			bestScore = mng.fillGapMinScore;
		} else {
			bound.setLeft(bestStart);
			for (int r = 0; r < leftCleavage.size(); r++) {
				double shifts[] = new double[1];
				shifts[0] = leftCleavage.get(r).getMass();
				ArrayList<EnumEnzyme> enzymeList = leftCleavage.get(r)
						.getEnzymeList();
				for (int i = 0; i < enzymeList.size(); i++) {
					EnumEnzyme enzyme = enzymeList.get(i);
					int id = enzyme.getId();
					for (int j = 0; j < enzymeBuSeqs[id].length; j++) {
						enzymeBuSeqs[id][j]
								.findLimitedBestShift(
										convertPeakList(curPeakList,
												curPeakList.size()), shifts,
										mng.peakTolerance);
						if (enzymeBuSeqs[id][j].getBestScore() >= bestScore) {
							double masses[] = enzymeBuSeqs[id][j]
									.getShiftedPrms();
							enzymeBuSeqs[id][j].output();
							compareRef(masses, masses.length);
							for (int k = 0; k < masses.length; k++) {
								if (masses[k] > gap.getLeft()
										&& masses[k] < gap.getRight()) {
									AlignPeak alignPeak = new AlignPeak(
											masses[k], enzymeBuSeqs[id][j]);
									alignPeak.setLeftAnchor(true);
									leftNewPeakList.add(alignPeak);
								}
							}
						}
					}
				}
			}
		}

		ArrayList<Cleavage> rightCleavage = getRightCleavages(
				rightAnchorPeakList, bound.getRight());
		System.out.println("right cleavage mass size " + rightCleavage.size());
		if (rightCleavage.size() > 0) {
			stop = false;
			bound.setRight(rightCleavage.get(0).getMass());
		}
		ArrayList<AlignPeak> rightNewPeakList = new ArrayList<AlignPeak>();
		bestScore = -1;
		double bestEnd = Double.MAX_VALUE;

		for (int r = 0; r < rightCleavage.size(); r++) {
			logger.info("right cleavage mass " + rightCleavage.get(r).getMass());
			double shifts[] = new double[1];
			ArrayList<EnumEnzyme> enzymeList = rightCleavage.get(r)
					.getEnzymeList();
			for (int i = 0; i < enzymeList.size(); i++) {
				EnumEnzyme enzyme = enzymeList.get(i);
				int id = enzyme.getId();
				for (int j = 0; j < enzymeBuSeqs[id].length; j++) {
					shifts[0] = rightCleavage.get(r).getMass()
							- enzymeBuSeqs[id][j].getBpSpec().getResSeq()
									.getResMassSum();
					enzymeBuSeqs[id][j].findLimitedBestShift(
							convertPeakList(curPeakList, curPeakList.size()),
							shifts, mng.peakTolerance);
					if (enzymeBuSeqs[id][j].getBestScore() > mng.fillGapMinScore) {
						peptideCount++;
						enzymeBuSeqs[id][j].output();
						System.out.println("Best shift "
								+ enzymeBuSeqs[id][j].getBestShift()
								+ " best score "
								+ enzymeBuSeqs[id][j].getBestScore());
					}
					if (enzymeBuSeqs[id][j].getBestScore() > bestScore) {
						bestScore = enzymeBuSeqs[id][j].getBestScore();
						bestEnd = rightCleavage.get(r).getMass();
					}
				}
			}
		}
		logger.debug("best right score " + bestScore);
		if (bestScore < mng.fillGapMinScore) {
			bestScore = mng.fillGapMinScore;
		} else {
			bound.setRight(bestEnd);
			for (int r = 0; r < rightCleavage.size(); r++) {
				double shifts[] = new double[1];
				ArrayList<EnumEnzyme> enzymeList = rightCleavage.get(r)
						.getEnzymeList();
				for (int i = 0; i < enzymeList.size(); i++) {
					EnumEnzyme enzyme = enzymeList.get(i);
					int id = enzyme.getId();
					for (int j = 0; j < enzymeBuSeqs[id].length; j++) {
						shifts[0] = rightCleavage.get(r).getMass()
								- enzymeBuSeqs[id][j].getBpSpec().getResSeq()
										.getResMassSum();
						enzymeBuSeqs[id][j]
								.findLimitedBestShift(
										convertPeakList(curPeakList,
												curPeakList.size()), shifts,
										mng.peakTolerance);
						if (enzymeBuSeqs[id][j].getBestScore() >= bestScore) {
							double masses[] = enzymeBuSeqs[id][j]
									.getShiftedPrms();
							enzymeBuSeqs[id][j].output();
							compareRef(masses, masses.length);
							for (int k = 0; k < masses.length; k++) {
								if (masses[k] > gap.getLeft()
										&& masses[k] < gap.getRight()) {
									AlignPeak alignPeak = new AlignPeak(
											masses[k], enzymeBuSeqs[id][j]);
									alignPeak.setRightAnchor(true);
									rightNewPeakList.add(alignPeak);
								}
							}
						}
					}
				}
			}
		}

		curPeakList.addAll(leftNewPeakList);
		curPeakList.addAll(rightNewPeakList);
		ArrayList<AlignPeak> newPeakList = removeSimilarPeak(curPeakList,
				new TdBuPeakComparator());
		Collections.sort(newPeakList, new PositionComparator());
		logger.info("EM completed");
		outputPeakList(newPeakList);
		logger.info("Keep local top ");
		newPeakList = keepLocalTop(newPeakList, bound.getLeft(),
				bound.getRight());
		outputPeakList(newPeakList);
		curPeakList.clear();
		curPeakList.addAll(newPeakList);
		compareRef(curPeakList, curPeakList.size());
		if (bound.getLeft() >= bound.getRight()) {
			stop = true;
		}
		return stop;
	}

	public void compareRef(double oriMasses[], int peakNum) {
		int count = 0;
		double masses[] = new double[peakNum];
		for (int i = 0; i < peakNum; i++) {
			masses[i] = oriMasses[i];
		}
		Arrays.sort(masses);
		String seq = refSeqSpec.getResSeq().getAcidString();
		double refPrms[] = refSeqSpec.getExtBMasses();
		boolean isMatches[] = new boolean[refPrms.length];
		Arrays.fill(isMatches, false);
		for (int i = 0; i < refPrms.length; i++) {
			boolean found = false;
			for (int j = 0; j < peakNum; j++) {
				double tolerance = mng.peakTolerance.compRelaxErrorTole(
						refPrms[i], precMassMinusWater);
				if (Math.abs(masses[j] - refPrms[i]) <= tolerance) {
					found = true;
					count++;
					break;
				}
			}
			isMatches[i] = found;
		}
		System.out.println("Total peak " + masses.length + " matched peak: "
				+ count);
		String anno = "";
		for (int i = 0; i < seq.length(); i++) {
			anno = anno + seq.charAt(i);
			if (isMatches[i + 1]) {
				anno = anno + "|";
			} else {
				anno = anno + " ";
			}
		}
		System.out.println("Anno " + anno);
	}

}
