package edu.iupui.proteomics.denovo3;

import java.util.Comparator;

public class LocalTopComparator implements Comparator<AlignPeak> {
	public int compare(AlignPeak p1, AlignPeak p2) {
		if (p1.getLocalTopScore() > p2.getLocalTopScore()) {
			return -1;
		} else if (p1.getLocalTopScore() < p2.getLocalTopScore()) {
			return 1;
		} else {
			return 0;
		}
	}
}
