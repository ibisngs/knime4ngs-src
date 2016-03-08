package de.helmholtz_muenchen.ibis.ngs.vcfutils;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "VCFutils" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class VCFutilsNodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the VCFutils node.
     */
    protected VCFutilsNodeDialog() {
    	
    	final SettingsModelString vcf = new SettingsModelString(VCFutilsNodeModel.CFGKEY_VCFFILE,null);
    	final SettingsModelString refvcf = new SettingsModelString(VCFutilsNodeModel.CFGKEY_REFVCFFILE,null);
    	final SettingsModelString utility = new SettingsModelString(VCFutilsNodeModel.CFGKEY_UTILITY,"");
    	final SettingsModelBoolean qstatsuseref = new SettingsModelBoolean(VCFutilsNodeModel.CFGKEY_QSTATSUSEREF, false);
    	final SettingsModelString snp = new SettingsModelString(VCFutilsNodeModel.CFGKEY_SNPFILE,null);
    	final SettingsModelString hapmap = new SettingsModelString(VCFutilsNodeModel.CFGKEY_HAPMAPFILE,null);
    	final SettingsModelIntegerBounded minrms = new SettingsModelIntegerBounded(VCFutilsNodeModel.CFGKEY_MINRMS, VCFutilsNodeModel.DEFAULT_MINRMS, 0, Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded minreaddepth = new SettingsModelIntegerBounded(VCFutilsNodeModel.CFGKEY_MINREADDEPTH, VCFutilsNodeModel.DEFAULT_MINREADDEPTH, 0, Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded maxreaddepth = new SettingsModelIntegerBounded(VCFutilsNodeModel.CFGKEY_MAXREADDEPTH, VCFutilsNodeModel.DEFAULT_MAXREADDEPTH, 0, Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded minaltbase = new SettingsModelIntegerBounded(VCFutilsNodeModel.CFGKEY_MINALTBASE, VCFutilsNodeModel.DEFAULT_MINALTBASE, 0, Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded gapfilter = new SettingsModelIntegerBounded(VCFutilsNodeModel.CFGKEY_GAPFILTER, VCFutilsNodeModel.DEFAULT_GAPFILTER, 0, Integer.MAX_VALUE);
    	final SettingsModelIntegerBounded adjacentgaps = new SettingsModelIntegerBounded(VCFutilsNodeModel.CFGKEY_ADJACENTGAPS, VCFutilsNodeModel.DEFAULT_ADJACENTGAPS, 0, Integer.MAX_VALUE);
    	final SettingsModelDoubleBounded strandpval = new SettingsModelDoubleBounded(VCFutilsNodeModel.CFGKEY_STRANDPVAL, VCFutilsNodeModel.DEFAULT_STRANDPVAL, 0, Double.MAX_VALUE);
    	final SettingsModelIntegerBounded baseqpval = new SettingsModelIntegerBounded(VCFutilsNodeModel.CFGKEY_BASEQPVAL, VCFutilsNodeModel.DEFAULT_BASEQPVAL, 0, Integer.MAX_VALUE);
    	final SettingsModelDoubleBounded mapqpval = new SettingsModelDoubleBounded(VCFutilsNodeModel.CFGKEY_MAPQPVAL, VCFutilsNodeModel.DEFAULT_MAPQPVAL, 0, Double.MAX_VALUE);
    	final SettingsModelDoubleBounded enddistpval = new SettingsModelDoubleBounded(VCFutilsNodeModel.CFGKEY_ENDDISTPVAL, VCFutilsNodeModel.DEFAULT_ENDDISTPVAL, 0, Double.MAX_VALUE);
    	final SettingsModelDoubleBounded hwepval = new SettingsModelDoubleBounded(VCFutilsNodeModel.CFGKEY_HWEPVAL, VCFutilsNodeModel.DEFAULT_HWEPVAL, 0, Double.MAX_VALUE);
    	final SettingsModelBoolean printfiltered = new SettingsModelBoolean(VCFutilsNodeModel.CFGKEY_PRINTFILTERED, false);
    	final SettingsModelIntegerBounded indelfilteringwindow = new SettingsModelIntegerBounded(VCFutilsNodeModel.CFGKEY_INDELFW, VCFutilsNodeModel.DEFAULT_INDELFW, 0, Integer.MAX_VALUE);

    	
    	qstatsuseref.setEnabled(false);
    	refvcf.setEnabled(false);
    	snp.setEnabled(false);
    	hapmap.setEnabled(false);
    	minrms.setEnabled(false);
    	minreaddepth.setEnabled(false);
    	maxreaddepth.setEnabled(false);
    	minaltbase.setEnabled(false);
    	gapfilter.setEnabled(false);
    	adjacentgaps.setEnabled(false);
    	strandpval.setEnabled(false);
    	baseqpval.setEnabled(false);
    	mapqpval.setEnabled(false);
    	enddistpval.setEnabled(false);
    	hwepval.setEnabled(false);
    	printfiltered.setEnabled(false);
    	indelfilteringwindow.setEnabled(false);
    	
    	createNewGroup("VCFutils");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(VCFutilsNodeModel.CFGKEY_VCF,null), "his_vcf_id", 0, ""));
    	addDialogComponent(new DialogComponentStringSelection(utility,"Select Utility", "subsam", "listsam", "fillac", "qstats", "hapmap2vcf", "ucscsnp2vcf", "varFilter", "vcf2fq"));
    	createNewGroup("Choose VCF file.");
    	addDialogComponent(new DialogComponentFileChooser(vcf, "his0_vcf_id", 0, "vcf","VCF"));
    	createNewGroup("Reference VCF file for qstats");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(qstatsuseref, "Use reference VCF file:"));
    	addDialogComponent(new DialogComponentFileChooser(refvcf, "his1_vcf_id", 0, "vcf","VCF"));
    	setHorizontalPlacement(false);
    	createNewGroup("SNP file for hapmap2vcf and ucscsnp2vcf");
    	addDialogComponent(new DialogComponentFileChooser(snp, "his2_vcf_id", 0, "snp","SNP"));
    	createNewGroup("HAPMAP file for hapmap2vcf");
    	addDialogComponent(new DialogComponentFileChooser(hapmap, "his3_vcf_id", 0, "hapmap","HAPMAP"));
    	createNewTab("Parameter");
    	createNewGroup("Parameter");
    	addDialogComponent(new DialogComponentNumber(minrms,"Minimum RMS mapping quality for SNPs:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(minreaddepth,"Minimum read depth:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(maxreaddepth,"Maximum read depth:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(minaltbase,"Minimum number of alternate bases:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(gapfilter,"SNP within x bp around a gap to be filtered:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(adjacentgaps,"Window size for filtering adjacent gaps:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(strandpval,"Min P-value for strand bias:", /*step*/ 0.0001));
    	addDialogComponent(new DialogComponentNumber(baseqpval,"Min P-value for baseQ bias [1.0E-]:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(mapqpval,"Min P-value for mapQ bias:", /*step*/ 0.0001));
    	addDialogComponent(new DialogComponentNumber(enddistpval,"Min P-value for end distance bias:", /*step*/ 0.0001));
    	addDialogComponent(new DialogComponentNumber(hwepval,"Min P-value for HWE (plus F<0):", /*step*/ 0.0001));
    	addDialogComponent(new DialogComponentBoolean(printfiltered, "Print filtered variants"));
    	addDialogComponent(new DialogComponentNumber(indelfilteringwindow,"INDEL filtering window:", /*step*/ 1));

    	
    	
    	
    	utility.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				String ut = utility.getStringValue();
				if(!VCFutilsNodeModel.optionalPort&&(ut.equals("subsam") || ut.equals("listsam") || ut.equals("fillac") || ut.equals("qstats") || ut.equals("varFilter") || ut.equals("vcf2fq"))) {
					vcf.setEnabled(true);
				} else {
					vcf.setEnabled(false);
				}
				qstatsuseref.setEnabled(false);
				refvcf.setEnabled(false);
				if(ut.equals("qstats")) {
					qstatsuseref.setEnabled(true);
					if(qstatsuseref.getBooleanValue()) {
						refvcf.setEnabled(true);
					}
				}
				if(ut.equals("hapmap2vcf")) {
					snp.setEnabled(true);
					hapmap.setEnabled(true);
				} else {
					snp.setEnabled(false);
					hapmap.setEnabled(false);
				}
				if(ut.equals("ucscsnp2vcf")) {
					snp.setEnabled(true);
				}
				if(ut.equals("varFilter")) {
					minrms.setEnabled(true);
			    	minreaddepth.setEnabled(true);
			    	maxreaddepth.setEnabled(true);
			    	minaltbase.setEnabled(true);
			    	gapfilter.setEnabled(true);
			    	adjacentgaps.setEnabled(true);
			    	strandpval.setEnabled(true);
			    	baseqpval.setEnabled(true);
			    	mapqpval.setEnabled(true);
			    	enddistpval.setEnabled(true);
			    	hwepval.setEnabled(true);
			    	printfiltered.setEnabled(true);
			    	minreaddepth.setIntValue(2);
			    	maxreaddepth.setIntValue(10000000);
				} else {
					minrms.setEnabled(false);
			    	minreaddepth.setEnabled(false);
			    	maxreaddepth.setEnabled(false);
			    	minaltbase.setEnabled(false);
			    	gapfilter.setEnabled(false);
			    	adjacentgaps.setEnabled(false);
			    	strandpval.setEnabled(false);
			    	baseqpval.setEnabled(false);
			    	mapqpval.setEnabled(false);
			    	enddistpval.setEnabled(false);
			    	hwepval.setEnabled(false);
			    	printfiltered.setEnabled(false);
				}
				if(ut.equals("vcf2fq")) {
					minreaddepth.setEnabled(true);
			    	maxreaddepth.setEnabled(true);
					minrms.setEnabled(true);
			    	indelfilteringwindow.setEnabled(true);
			    	minreaddepth.setIntValue(3);
			    	maxreaddepth.setIntValue(100000);
				}else{
					indelfilteringwindow.setEnabled(false);
				}
			}
		});
    	
    	qstatsuseref.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				refvcf.setEnabled(qstatsuseref.getBooleanValue());
			}
		});
    	
    	
    }
}

