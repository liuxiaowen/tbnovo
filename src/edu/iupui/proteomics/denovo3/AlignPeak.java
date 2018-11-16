package edu.iupui.proteomics.denovo3;

import java.util.ArrayList;

import edu.iupui.proteomics.spec.peak.Peak;

public class AlignPeak implements Peak{
	private double position;
	private double intensity;
	private ArrayList<TdDenovoPeak> tdPeakList;
	private ArrayList<BottomUpSeq> buSeqList;
	
	private AlignPeak prevPeak = null;
	private double curScore = 0;
	private int length = 0;
	
	/* isLeftMatch, isExactMatch and isRightMatch are used for annotation */
	private boolean isLeftMatch = false;
	private boolean isExactMatch = false;
	private boolean isRightMatch = false;

	/* isLeftAnchor and isRightAnchor are used for gap filling */
	private boolean isLeftAnchor = false;
	private boolean isRightAnchor = false;
	
	/* the confidence score of top-down spectra */
	private double tdConf = 0.5;
//	public static double LOW_TD_CONF = 0.5;
//	public static double NORMAL_TD_CONF = 1.0;
			

    public AlignPeak(TdDenovoPeak tdPeak) {
    	tdPeakList = new ArrayList<TdDenovoPeak> ();
    	tdPeakList.add(tdPeak);
    	this.position = tdPeak.getPosition();
    	this.intensity = tdPeak.getIntensity();
    	buSeqList = new ArrayList<BottomUpSeq>();
    }
    
    public AlignPeak(AlignPeak peak) {
    	position = peak.getPosition();
    	intensity = peak.getIntensity();
    	tdPeakList = peak.getTdPeakList();
    	buSeqList = new ArrayList<BottomUpSeq>();
    }
    
    public AlignPeak(double pos, BottomUpSeq seq) {
    	this.position = pos;
    	tdPeakList = new ArrayList<TdDenovoPeak>();
    	buSeqList = new ArrayList<BottomUpSeq>();
    	buSeqList.add(seq);
    }
    
    public AlignPeak(double pos) {
    	this.position = pos;
    	tdPeakList = new ArrayList<TdDenovoPeak>();
    	buSeqList = new ArrayList<BottomUpSeq>();
    }
    
    public double getPosition() {
        return position;
    }
       
    public void setPosition(double position) {
        this.position = position;
    }

	public ArrayList<TdDenovoPeak> getTdPeakList() {
		return tdPeakList;
	}

	public int getTdPeakLen() {
		return tdPeakList.size();
	}
	
	public double getIntensity() {
		return intensity;
	}

	public ArrayList<BottomUpSeq> getSeqList() {
		return buSeqList;
	}
	
	public int getBuSeqLen() {
		return buSeqList.size();
	}
	
	public void mergePeak(AlignPeak peak) {
		ArrayList<TdDenovoPeak> newTdPeakList = peak.getTdPeakList();
		for (int i = 0; i < newTdPeakList.size(); i++) {
			TdDenovoPeak newPeak = newTdPeakList.get(i);
			if (!tdPeakList.contains(newPeak)) {
				tdPeakList.add(newPeak);
			}
		}
		ArrayList<BottomUpSeq> newSeqList = peak.getSeqList();
		for (int i = 0; i < newSeqList.size(); i++) {
			BottomUpSeq newSeq = newSeqList.get(i);
			if (!buSeqList.contains(newSeq)) {
				buSeqList.add(newSeq);
			}
		}
	}

	public void setIntensity(double intensity) {
		this.intensity = intensity;
	}

	public double getAlignScore() {
		double tdScore = 0;
		if (tdPeakList.size() != 0) {
			tdScore  = Math.log(tdPeakList.size())/Math.log(2) + 1;
		}
		double buScore = 0;
		for (int i = 0; i < buSeqList.size(); i++) {
			buScore = buScore + Math.log(buSeqList.get(i).getWeight())/Math.log(2) + 1;
		}
		return tdScore + buScore;
	}
	
	public double getLocalTopScore() {
		double tdScore = 0;
		if (tdPeakList.size() != 0) {
			tdScore  = tdConf;
		}
		double buScore = 0;
		for (int i = 0; i < buSeqList.size(); i++) {
			buScore = buScore + Math.log(buSeqList.get(i).getWeight())/Math.log(2) + 1;
		}
		return tdScore + buScore;
	}
	
	public void removeBuSeqs(ArrayList<BottomUpSeq> seqList) {
		BottomUpSeq seqs[] = buSeqList.toArray(new BottomUpSeq[0]);
		for (int i =0; i < seqs.length; i++) {
			if (seqList.contains(seqs[i])) {
				buSeqList.remove(seqs[i]);
			}
		}
	}
	
	public void removeAllBuSeqs() {
		buSeqList.clear();
	}

	public AlignPeak getPrevPeak() {
		return prevPeak;
	}

	public void setPrevPeak(AlignPeak prevPeak) {
		this.prevPeak = prevPeak;
	}

	public double getCurScore() {
		return curScore;
	}

	public void setCurScore(double curScore) {
		this.curScore = curScore;
	}


	public boolean isLeftMatch() {
		return isLeftMatch;
	}

	public void setLeftMatch(boolean isLeftMatch) {
		this.isLeftMatch = isLeftMatch;
	}

	public boolean isExactMatch() {
		return isExactMatch;
	}

	public void setExactMatch(boolean isExactMatch) {
		this.isExactMatch = isExactMatch;
	}

	public boolean isRightMatch() {
		return isRightMatch;
	}

	public void setRightMatch(boolean isRightMatch) {
		this.isRightMatch = isRightMatch;
	}
	
	public boolean isLeftAnchor() {
		return isLeftAnchor;
	}

	public void setLeftAnchor(boolean isLeftAnchor) {
		this.isLeftAnchor = isLeftAnchor;
	}

	public boolean isRightAnchor() {
		return isRightAnchor;
	}

	public void setRightAnchor(boolean isRightAnchor) {
		this.isRightAnchor = isRightAnchor;
	}

	public double getTdConf() {
		return tdConf;
	}

	public void setTdConf(double tdConf) {
		this.tdConf = tdConf;
	}

	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}
}
