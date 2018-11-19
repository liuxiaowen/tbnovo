package edu.iupui.proteomics.denovo.reader;

import java.util.ArrayList;

import edu.iupui.proteomics.base.seq.ResSeq;
import edu.iupui.proteomics.base.theosp.BpSpec;

public class DenovoBpSpecReader {
	/** initialize sequence list */
	public static BpSpec[] readDb(DenovoSeqReader reader) throws Exception {
		ArrayList<BpSpec> bpSpecList = new ArrayList<BpSpec>();
		ResSeq resSeq;
		int id = 0;
		while ((resSeq = reader.getNextResSeq()) != null) {
			resSeq.setId(id);
			id++;
			bpSpecList.add(new BpSpec(resSeq));
		}
		BpSpec bpSpecs[] = bpSpecList.toArray(new BpSpec[0]);
		bpSpecList = null;
		reader = null;
		return bpSpecs;
	}

}
