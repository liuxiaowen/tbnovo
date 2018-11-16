package edu.iupui.proteomics.denovo3;

import java.util.Comparator;


public class NormScoreComparator implements Comparator<BottomUpSeq> {
	public int compare(BottomUpSeq s1, BottomUpSeq s2) {
		if (s1.getNormalizedScore() > s2.getNormalizedScore()) {
			return -1;
		} else if (s1.getNormalizedScore() < s2.getNormalizedScore()) {
			return 1;
		} else {
			return 0;
		}
	}
}
