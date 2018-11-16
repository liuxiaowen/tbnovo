package edu.iupui.proteomics.denovo3;

import java.io.File;

import edu.iupui.proteomics.base.residue.ResList;
import edu.iupui.proteomics.base.residue.ResListFactory;
import edu.iupui.proteomics.base.theosp.BpSpec;
import edu.iupui.proteomics.spec.deconvsp.DeconvPeak;
import edu.iupui.proteomics.spec.deconvsp.reader.MsAlignReader;
import edu.iupui.proteomics.spec.sp.Ms;

import edu.iupui.proteomics.denovo.reader.DenovoBpSpecReader;
import edu.iupui.proteomics.denovo.reader.DenovoSeqReader;

public class DenovoProcessor {
    private DenovoMng mng;
    
    public DenovoProcessor(DenovoMng mng) {
    	this.mng = mng;
    }
    
	public void process() throws Exception{
		/* read top-down spectrum */ 
        File tdSpFile = new File(mng.tdSpFileName);
        int tdSpNum = MsAlignReader.countSpNum(tdSpFile);
        System.out.println("Td spectra number " + tdSpNum);
        MsAlignReader tdSpReader = new MsAlignReader(tdSpFile);

        @SuppressWarnings("unchecked")
		Ms<DeconvPeak> tdSps[] = new Ms[tdSpNum];
        @SuppressWarnings("unchecked")
		Ms<TdDenovoPeak> tdPrmSps[] = new Ms[tdSpNum];
        for (int i = 0; i < tdSpNum; i++) {
		    Ms<DeconvPeak>[] sps = tdSpReader.getNextMses();
            tdSps[i] = sps[0];
            tdPrmSps[i] = TdDenovoMsFactory.getMsTwo(tdSps[i], 0, mng.minMass,
				mng.peakTolerance);
        }
        tdSpReader.close();
        
        /* read bottom-up spectrum */
        File denovoFile = new File(mng.buSpFileName);
        ResList buResList = ResListFactory.getSystemInstance(mng.buResFileName);
        BottomUpSeq buSeqs[] = BottomUpSeqReader.readDb(new DenovoSeqReader(denovoFile, buResList));
 
        File refFile = new File(mng.refSeqFileName);
        ResList resList = ResListFactory.getSystemInstance(mng.refResFileName);
        BpSpec refSeqs[] = DenovoBpSpecReader.readDb(new DenovoSeqReader(refFile, resList));
    
        DenovoAlign align = new DenovoAlign(mng, tdPrmSps, buSeqs, refSeqs[0]);
        align.process();
	}
}
