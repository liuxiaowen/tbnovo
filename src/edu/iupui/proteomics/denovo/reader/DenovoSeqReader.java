package edu.iupui.proteomics.denovo.reader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;

import edu.iupui.proteomics.base.residue.Acid;
import edu.iupui.proteomics.base.residue.AcidList;
import edu.iupui.proteomics.base.residue.Res;
import edu.iupui.proteomics.base.residue.ResArrayUtil;
import edu.iupui.proteomics.base.residue.ResList;
import edu.iupui.proteomics.base.seq.ResSeq;

public class DenovoSeqReader {
	BufferedReader input;
	public AcidList acidList;
	public ResList resList;

	/**
	 * Constructs an instance with a File.
	 */
	public DenovoSeqReader(File file, ResList resList) throws Exception {
		input = new BufferedReader(new InputStreamReader(new FileInputStream(
				file), "UTF-8"));
		this.acidList = AcidList.getCompleteAcidList();
		this.resList = resList;
	}

	/**
	 * Read FASTA file and return next protein name and sequence. result[0] is
	 * protein name and result[1] is sequence.
	 */
	public String[] getNextSeq() throws Exception {
		if (input == null) {
			return null;
		}
		String line = input.readLine();
		if (line == null) {
			input = null;
			return null;
		}
		String results[] = new String[2];
		String words[] = line.split(",");
		System.out.println(line + " word 0 " + words[0]);
		results[0] = words[1];
		results[1] = words[3];
		return results;
	}
	
	public String[] getNextLine() throws Exception {
		if (input == null) {
			return null;
		}
		String line = input.readLine();
		if (line == null) {
			input = null;
			return null;
		}
		String words[] = line.split(",");
		return words;
	}

	/**
	 * Read FASTA file and return next protein as an ResSeq.
	 */
	public ResSeq getNextResSeq() throws Exception {
		String[] seqInfo = getNextSeq();
		if (seqInfo == null) {
			return null;
		}
		String name = seqInfo[0];
		String seq = seqInfo[1];
		Acid acids[] = acidList.convert(seq);
		Res residues[] = ResArrayUtil.getResArrayByAcid(resList, acids);
		return new ResSeq(name, residues);
	}

}
