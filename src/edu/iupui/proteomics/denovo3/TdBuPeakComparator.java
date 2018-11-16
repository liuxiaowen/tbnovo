package edu.iupui.proteomics.denovo3;

import java.util.Comparator;

public class TdBuPeakComparator implements Comparator<AlignPeak> {
	public int compare(AlignPeak p1, AlignPeak p2) {
		if (p1.getTdPeakLen() > 0 && p2.getTdPeakLen() == 0) {
			return -1;
		} else if (p1.getTdPeakLen() == 0 && p2.getTdPeakLen() > 0) {
			return 1;
		} else {
			if (p1.getBuSeqLen() > p2.getBuSeqLen()) {
				return -1;
			} else if (p1.getBuSeqLen() < p2.getBuSeqLen()) {
				return 1;
			} else {
				if (p1.getTdPeakLen() > p2.getTdPeakLen()) {
					return -1;
				} else if (p1.getTdPeakLen() < p2.getTdPeakLen()) {
					return 1;
				} else {

					if (p1.getIntensity() > p2.getIntensity()) {
						return -1;
					} else if (p1.getIntensity() < p2.getIntensity()) {
						return 1;
					} else {
						return 0;
					}
				}
			}
		}
	}
}
