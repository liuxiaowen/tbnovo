package edu.iupui.proteomics.denovo3;

import java.util.ArrayList;
import java.util.Collections;
import java.util.SortedMap;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import edu.iupui.proteomics.base.residue.Res;
import edu.iupui.proteomics.base.residue.ResList;
import edu.iupui.proteomics.base.residue.ResListFactory;
import edu.iupui.proteomics.spec.peak.PositionComparator;

public class FindPath {

	public static Logger logger = Logger.getLogger(FindPath.class);

	public ResList resList;
	public SortedMap<Double, Res[]> massMap;
	public ArrayList<AlignPeak> peakList;
	private double precMassMinusWater;
	private DenovoMng mng;

	private ArrayList<AlignPeak> resultPeakList;

	public FindPath(ArrayList<AlignPeak> peakList, double precMassMinusWater,
			DenovoMng mng) throws Exception {
		resList = ResListFactory.getSystemInstance(mng.buResFileName);
		logger.debug("ResList length " + resList.size());
		this.peakList = peakList;
		this.precMassMinusWater = precMassMinusWater;
		this.mng = mng;

		massMap = new TreeMap<Double, Res[]>();
		for (int i = 0; i < 20; i++) {
			double mass = resList.get(i).getMass();
			Res residues[] = new Res[1];
			residues[0] = resList.get(i);
			if (!massMap.containsKey(mass)) {
				massMap.put(mass, residues);
			}
		}
		for (int i = 0; i < 20; i++) {
			for (int j = i + 1; j < 20; j++) {
				double mass = resList.get(i).getMass()
						+ resList.get(j).getMass();
				Res residues[] = new Res[2];
				residues[0] = resList.get(i);
				residues[1] = resList.get(j);
				if (!massMap.containsKey(mass)) {
					massMap.put(mass, residues);
				}
			}
		}

		// for (int i = 0; i < 20; i++) {
		// for (int j = i + 1; j < 20; j++) {
		// for (int k = j + 1; k < 20; k++) {
		// double mass = resList.get(i).getMass()
		// + resList.get(j).getMass()
		// + resList.get(k).getMass();
		// Res residues[] = new Res[3];
		// residues[0] = resList.get(i);
		// residues[1] = resList.get(j);
		// residues[2] = resList.get(k);
		// if (!massMap.containsKey(mass)) {
		// massMap.put(mass, residues);
		// }
		// }
		// }
		// }
	}

	public Res[] process() {
		dp();
		Res residues[] = backTrack();
		return residues;
	}

	private void dp() {
		peakList.add(new AlignPeak(0));
		peakList.add(new AlignPeak(precMassMinusWater));
		Collections.sort(peakList, new PositionComparator());
		for (int i = 1; i < peakList.size(); i++) {
			AlignPeak curPeak = peakList.get(i);
			AlignPeak bestPrevPeak = null;
			double bestScore = Double.NEGATIVE_INFINITY;
			for (int j = i - 1; j >= 0; j--) {
				AlignPeak prevPeak = peakList.get(j);
				double diff = curPeak.getPosition() - prevPeak.getPosition();
				if (diff > 600) {
					break;
				}
				Double key = find(diff);
				if (key != null) {
					Res residues[] = massMap.get(key);
					double newScore = prevPeak.getCurScore() - residues.length
							* 2.5;
					if (newScore > bestScore) {
						bestPrevPeak = prevPeak;
						bestScore = newScore;
					}
				}
			}

			if (bestPrevPeak != null) {
				curPeak.setCurScore(bestScore + curPeak.getAlignScore());
				curPeak.setLength(bestPrevPeak.getLength() + 1);
				System.out.println("mass " + curPeak.getPosition()
						+ " pre mass " + bestPrevPeak.getPosition() + " score "
						+ curPeak.getAlignScore() + " total score "
						+ curPeak.getCurScore());
			} else {
				int bestLength = -1;
				for (int j = 0; j < i; j++) {
					AlignPeak prevPeak = peakList.get(j);
					int newLength = prevPeak.getLength();
					if (newLength > bestLength) {
						bestPrevPeak = prevPeak;
						bestLength = newLength;
					}

				}
				curPeak.setCurScore(0);
				curPeak.setLength(bestPrevPeak.getLength());
				System.out.println("mass " + curPeak.getPosition()
						+ " pre mass null ");
			}
			curPeak.setPrevPeak(bestPrevPeak);
		}
	}

	private Res[] backTrack() {
		ArrayList<Res> residueList = new ArrayList<Res>();
		AlignPeak peak = peakList.get(peakList.size() - 1);
		resultPeakList = new ArrayList<AlignPeak>();
		while (peak.getPrevPeak() != null) {
			logger.debug("Mass " + peak.getPosition() + " score "
					+ peak.getAlignScore() + " bu number " + peak.getBuSeqLen());
			double diff = peak.getPosition() - peak.getPrevPeak().getPosition();
			Double key = find(diff);
			if (key != null) {
				Res residues[] = massMap.get(key);
				for (int i = residues.length - 1; i >= 0; i--) {
					residueList.add(residues[i]);
					System.out.print(residues[i].getAcid().getOneLetter());
				}
				System.out.println();
			}
			peak = peak.getPrevPeak();
			resultPeakList.add(peak);
		}
		Res results[] = new Res[residueList.size()];
		for (int i = 0; i < results.length; i++) {
			results[i] = residueList.get(results.length - i - 1);
		}
		return results;
	}

	private Double find(double diff) {
		double minErrorOne = 1000;
		double minErrorTwo = 1000;
		Double bestKeyOne = null;
		Double bestKeyTwo = null;
		for (Double key : massMap.keySet()) {
			if (Math.abs(diff - key) < minErrorOne
					&& massMap.get(key).length == 1) {
				minErrorOne = Math.abs(diff - key);
				bestKeyOne = key;
			}
		}

		for (Double key : massMap.keySet()) {
			if (Math.abs(diff - key) < minErrorTwo
					&& massMap.get(key).length == 2) {
				minErrorTwo = Math.abs(diff - key);
				bestKeyTwo = key;
			}
		}
		if (minErrorOne < mng.findPathMinError) {
			return bestKeyOne;
		}
		if (minErrorTwo < mng.findPathMinError) {
			return bestKeyTwo;
		}
		return null;
	}

	public ArrayList<AlignPeak> getResultPeakList() {
		return resultPeakList;
	}

}
