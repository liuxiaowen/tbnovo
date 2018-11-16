package edu.iupui.proteomics.denovo3;

import edu.iupui.proteomics.base.util.BioIo;

public class DenovoConsole {
	public static void main(String args[]) {
        try {
        	DenovoMng mng = new DenovoMng();
            mng.tdSpFileName = args[0];

            mng.buSpFileName = args[1];
            mng.buResFileName = args[2];

            mng.refSeqFileName = args[3];
            mng.refResFileName = args[4];
            
            mng.minBuScore = Double.parseDouble(args[5]);
            
            mng.fillGapMinScore = Double.parseDouble(args[6]);
            
            mng.seleSeqNum = Integer.parseInt(args[7]);
            
            mng.minSupport = Integer.parseInt(args[8]);
            
            mng.confBuSeqThresh = Integer.parseInt(args[9]);
            
            DenovoProcessor processor = new DenovoProcessor(mng);
            processor.process();
        } catch (Exception e) {
            System.out.println("\n\n" + BioIo.getStackTraceAsString(e));
            System.exit(1);
        } catch (Error e) {
            System.out.println("\n\n" + BioIo.getStackTraceAsString(e));
            System.exit(1);
        }

	}

}
