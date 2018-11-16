package edu.iupui.proteomics.denovo3;

import edu.iupui.proteomics.base.ion.EnumActivation;
import edu.iupui.proteomics.spec.extendsp.ExtendSpPara;
import edu.iupui.proteomics.spec.peak.PeakTolerance;
import edu.iupui.proteomics.spec.sp.SpPara;

public class DenovoMng {

	public String tdSpFileName;
	
	public String buSpFileName;
	public String buResFileName;
	
	public String refSeqFileName;
	public String refResFileName;
	
	
    public double minMass = 50f;
    /** error tolerance */
    private double ppo = 0.000015f;
    private double minTolerance = 0.025f;
    /** class for computing error tolerance for peaks */
    public PeakTolerance peakTolerance = new PeakTolerance(ppo, true, minTolerance);
    
    public double convertRatio = 274.335215;
    
    /** top down peak filtration */
    
    /** error tolerance for merging top-down prm peaks */
    public double peakMergeTolerance = 0.2;
    /** top down peaks with only one supporting peak is removed */
    public int minSupport = 3;
    
    /** minimum bottom-up top-down matching score for bottom-up spectra */
    public double confBuSeqThresh = 70;
    public double minBuScore = 6.5;
    
    public double gapLen = 300;
    public double flankLen = 3000;
    public double fillgGapTolerance = 0.15;
    public double fillGapMinScore = 4.5;
    public int keepLocalNum = 3;
    
    public double longGapLen = 500;
    public double longGapPrecTolerance = 0.50;
    
    public int seleSeqNum = 3;

    /// continue
    
    public double minError = 0.05;
    public int size2 = 91;
    public int size3 = 67;
    //public int size5 = 100;
    
    public double findPathMinError = 0.2;
    
    private int minPeakNum = 0;
    /** the set of offsets used to expand the monoisotopic mass list */
    private double extOffsets[] = { 0f};
    private double extendThresh = 5000;
    private ExtendSpPara extendSpPara = new ExtendSpPara(extendThresh,
            extOffsets);
    private EnumActivation activationType = EnumActivation.CID;
    public  SpPara spPara = new SpPara(minPeakNum, minMass, peakTolerance, extendSpPara, activationType);
    
    public int convertToInt(double value) {
    	return (int)Math.round(value * convertRatio);
    }

}
