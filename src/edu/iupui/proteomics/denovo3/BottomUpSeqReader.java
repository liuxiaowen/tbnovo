package edu.iupui.proteomics.denovo3;

import java.io.File;
import java.util.ArrayList;

import edu.iupui.proteomics.base.residue.Acid;
import edu.iupui.proteomics.base.residue.Res;
import edu.iupui.proteomics.base.residue.ResArrayUtil;
import edu.iupui.proteomics.base.seq.ResSeq;
import edu.iupui.proteomics.base.theosp.BpSpec;
import edu.iupui.proteomics.spec.deconvsp.DeconvPeak;
import edu.iupui.proteomics.spec.deconvsp.DeconvSpFactory;
import edu.iupui.proteomics.spec.deconvsp.reader.MsAlignReader;
import edu.iupui.proteomics.spec.prmsp.PrmPeak;
import edu.iupui.proteomics.spec.prmsp.PrmSpFactory;
import edu.iupui.proteomics.spec.rawsp.RawPeak;
import edu.iupui.proteomics.spec.rawsp.simplereader.MgfSimpleReader;
import edu.iupui.proteomics.spec.rawsp.simplereader.MsSimpleReader;
import edu.iupui.proteomics.spec.sp.Ms;
import edu.iupui.proteomics.spec.sp.SpPara;
import edu.iupui.proteomics.spec.spset.SpectrumSet;

import edu.iupui.proteomics.denovo.reader.DenovoSeqReader;
import edu.ucsd.proteomics.msdeconv.sp.reader.MsReader;
import edu.ucsd.proteomics.msdeconv.sp.reader.MzXmlSingleReader;

public class BottomUpSeqReader {
	public static BottomUpSeq[] readDb(DenovoSeqReader reader) throws Exception {
		ArrayList<BottomUpSeq> buSeqList = new ArrayList<BottomUpSeq>();
		String words[];
		int id = 0;
		while ((words = reader.getNextLine()) != null) {
			String name = words[1];
			String raw_seq = words[3];
			String parts[] = raw_seq.split(" ");
			double nGap = 0.0;
			double cGap = 0.0;
			String seq;
			if (parts.length == 1) {
				seq = raw_seq;
			}
			else {
				nGap = Double.parseDouble(parts[0]);
				cGap = Double.parseDouble(parts[2]);
				seq = parts[1];
			}
		
			Acid acids[] = reader.acidList.convert(seq);
			Res residues[] = ResArrayUtil.getResArrayByAcid(reader.resList,
					acids);
			ResSeq resSeq = new ResSeq(name, residues);
			resSeq.setId(id);
			id++;
			int scan = Integer.parseInt(words[2]);
			String fileName = words[0];
			double denovoScore = Double.parseDouble(words[5]);
			BottomUpSeq buSeq = new BottomUpSeq(new BpSpec(resSeq), nGap, cGap, fileName,
					scan, denovoScore);
			buSeqList.add(buSeq);
		}
		BottomUpSeq buSeqs[] = buSeqList.toArray(new BottomUpSeq[0]);
		return buSeqs;
	}

	/*
	public static void addMzXMLSpectra(BottomUpSeq seqs[], String fileName,
			SpPara spPara) throws Exception {
		String spFileName = fileName + ".mzXML";
		MsReader reader = (MsReader) (new MzXmlSingleReader(spFileName));
		Ms<RawPeak> sp;
		while ((sp = reader.getNextMs()) != null) {
			int scan = sp.getHeader().getFirstScanNum();
			for (int i = 0; i < seqs.length; i++) {
				if (seqs[i].getScan() == scan
						&& seqs[i].getFileName().equals(fileName)) {
					Ms<DeconvPeak> deconvSp = DeconvSpFactory.getDeconvMs(sp,
							spPara);
					double shift = seqs[i].getBpSpec().getResSeq().getSeqMass()
							- sp.getHeader().getPrecMonoMass();
					Ms<PrmPeak> prmSp = PrmSpFactory.getMsTwo(deconvSp, shift,
							spPara);
					seqs[i].setPrmSp(prmSp);
				}
			}

		}

	}
	*/
	
	public static void addSpectra(BottomUpSeq seqs[], String fileName,
			SpPara spPara) throws Exception {
		String spFileName = fileName + ".mgf";
		MsSimpleReader reader = new MgfSimpleReader(spFileName, true);
		Ms<RawPeak> sp;
		while ((sp = reader.getNextMs()) != null) {
			int scan = sp.getHeader().getFirstScanNum();
			for (int i = 0; i < seqs.length; i++) {
				if (seqs[i].getScan() == scan
						&& seqs[i].getFileName().equals(fileName)) {
					Ms<DeconvPeak> deconvSp = DeconvSpFactory.getDeconvMs(sp,
							spPara);
					double shift = seqs[i].getBpSpec().getResSeq().getSeqMass()
							- sp.getHeader().getPrecMonoMass();
					Ms<PrmPeak> prmSp = PrmSpFactory.getMsTwo(deconvSp, shift,
							spPara);
					seqs[i].setPrmSp(prmSp);
				}
			}

		}

	}

	public static void addMsdeconvSpectra(BottomUpSeq seqs[], String fileName,
			SpPara spPara) throws Exception {
		String spFileName = fileName + "_msdeconv.msalign";

		MsAlignReader spReader = new MsAlignReader(new File(spFileName));
		Ms<DeconvPeak>[] deconvSp;
		while ((deconvSp = spReader.getNextMses()) != null) {
			SpectrumSet spectrumSet = SpectrumSet.getSpectrumSet(deconvSp[0],
					0, spPara);
			if (spectrumSet == null) {
				System.out.println("null spectrumSet " + " file name " + fileName + " scan  " + deconvSp[0].getHeader().getFirstScanNum());
			}
			if (spectrumSet != null) {
				int scan = deconvSp[0].getHeader().getFirstScanNum();
				for (int i = 0; i < seqs.length; i++) {
					if (seqs[i].getScan() == scan
							&& seqs[i].getFileName().equals(fileName)) {
						Ms<PrmPeak> prmSp = spectrumSet.getSpTwo();
						seqs[i].setPrmSp(prmSp);
					}
				}

			}

		}
	}

}
