package edu.iupui.proteomics.denovo3;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.log4j.Logger;

import edu.iupui.proteomics.base.ion.EnumActivation;
import edu.iupui.proteomics.base.ion.EnumIonType;
import edu.iupui.proteomics.spec.deconvsp.DeconvPeak;
import edu.iupui.proteomics.spec.deconvsp.DeconvSpFactory;
import edu.iupui.proteomics.spec.peak.PeakTolerance;
import edu.iupui.proteomics.spec.peak.PositionComparator;
import edu.iupui.proteomics.spec.prmsp.EnumPrmPeakType;
import edu.iupui.proteomics.spec.sp.Ms;
import edu.iupui.proteomics.spec.sp.MsHeader;

public class TdDenovoMsFactory {
    private static Logger logger = Logger.getLogger(TdDenovoMsFactory.class);

    public static Ms<TdDenovoPeak> getMsTwo(Ms<DeconvPeak> deconvSp,
            double delta, double blankMass, PeakTolerance peakTolerance)
            throws Exception {
        logger.trace("start sp two generation");
        MsHeader header = DeconvSpFactory.getDeltaHeader(delta, deconvSp);
        DeconvPeak peaks[] = deconvSp.toArray(new DeconvPeak[0]);

        TdDenovoPeak denovoPeaks[] = getSpTwoDenovoPeak(peaks, header, blankMass);

        setTolerance(denovoPeaks, header.getPrecMonoMass(), peakTolerance);
        return new Ms<TdDenovoPeak>(denovoPeaks, header);
    }

    /*******************
     * private functions
     *******************/

    private static void addTwoMasses(ArrayList<TdDenovoPeak> list, DeconvPeak peak,
            double precMonoMass, EnumActivation activeType) throws Exception {
        /* generate extended vertex based on m */
        double origMass = peak.getMonoMass() - activeType.getNShift();
        TdDenovoPeak newPeak = new TdDenovoPeak(peak, origMass, EnumPrmPeakType.ORIGINAL, 1);
        list.add(newPeak);
        double reverseMass = precMonoMass
                - (peak.getMonoMass() - activeType.getCShift());
        newPeak = new TdDenovoPeak(peak, reverseMass, EnumPrmPeakType.REVERSED, 1);
        list.add(newPeak);

    }

    /**
     * generate type 2 PRM peak list: both original spectrum and reverse
     * spectrum. But +/-1 Da shift is not considered.
     */
    private static TdDenovoPeak[] getSpTwoDenovoPeak(DeconvPeak peaks[],
            MsHeader header, double blankMass) throws Exception {
        double precMonoMass = header.getPrecMonoMass();
        EnumActivation activeType = header.getActivationType();
        ArrayList<TdDenovoPeak> list = new ArrayList<TdDenovoPeak>();
        for (int i = 0; i < peaks.length; i++) {
            addTwoMasses(list, peaks[i], precMonoMass, activeType);
        }
        filterTdDenovoPeak(list, blankMass, precMonoMass);
        /* add 0 and precursor mass - 18 */
        DeconvPeak zeroPeak = new DeconvPeak(-1, 0, 0, 0);
        list.add(new TdDenovoPeak(zeroPeak, 0, EnumPrmPeakType.ORIGINAL, 1));
        DeconvPeak precPeak = new DeconvPeak(-1, precMonoMass
                - EnumIonType.PREC.getShift(), 0, 0);
        list.add(new TdDenovoPeak(precPeak, precMonoMass
                - EnumIonType.PREC.getShift(), EnumPrmPeakType.ORIGINAL, 1));
        /* sort by mass */
        TdDenovoPeak extPeaks[] = list.toArray(new TdDenovoPeak[0]);
        Arrays.sort(extPeaks, new PositionComparator());
        return extPeaks;
    }

    private static void setTolerance(TdDenovoPeak peaks[], double precMass,
            PeakTolerance peakTolerance) {
        peaks[0].setTolerance(peakTolerance.compStrictErrorTole(0));
        peaks[peaks.length - 1].setTolerance(peakTolerance
                .compStrictErrorTole(precMass));

        for (int i = 1; i < peaks.length - 1; i++) {
            /* if peak is m */
            peaks[i].setTolerance(peakTolerance
                    .compStrictErrorTole(peaks[i].getBasePeak().getMonoMass()));
        }
    }

    /** remove peaks near 0 and prec mass */
    private static void filterTdDenovoPeak(ArrayList<TdDenovoPeak> list,
            double blankMass, double precMass) {
        TdDenovoPeak peaks[] = list.toArray(new TdDenovoPeak[0]);
        for (int i = 0; i < peaks.length; i++) {
            double mass = peaks[i].getPosition();
            if (mass < blankMass || mass > precMass - blankMass) {
                list.remove(peaks[i]);
            }
        }
    }
}
