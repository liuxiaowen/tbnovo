package edu.iupui.proteomics.denovo3;

import java.util.Comparator;


public class DenovoScoreComparator implements Comparator<BottomUpSeq> {
	public int compare(BottomUpSeq s1, BottomUpSeq s2) {
		if (s1.getDenovoScore() > s2.getDenovoScore()) {
			return -1;
		} else if (s1.getDenovoScore() < s2.getDenovoScore()) {
			return 1;
		} else {
			return 0;
		}
	}
}
