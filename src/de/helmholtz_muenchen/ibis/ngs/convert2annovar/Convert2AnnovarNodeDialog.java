package de.helmholtz_muenchen.ibis.ngs.convert2annovar;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "Convert2Annovar" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class Convert2AnnovarNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the Convert2Annovar node.
     */
    protected Convert2AnnovarNodeDialog() {
    	
    	final SettingsModelString inFile = new SettingsModelString(Convert2AnnovarNodeModel.CFGKEY_INFILE,null);
    	final SettingsModelIntegerBounded snpqual = new SettingsModelIntegerBounded(Convert2AnnovarNodeModel.CFGKEY_SNPQUAL, Convert2AnnovarNodeModel.DEFAULT_SNPQUAL, 0, Integer.MAX_VALUE);
    	final SettingsModelDoubleBounded snppvalue = new SettingsModelDoubleBounded(Convert2AnnovarNodeModel.CFGKEY_SNPPVALUE, Convert2AnnovarNodeModel.DEFAULT_SNPPVALUE, 0, 1);
    	final SettingsModelIntegerBounded coverage = new SettingsModelIntegerBounded(Convert2AnnovarNodeModel.CFGKEY_COVERAGE, Convert2AnnovarNodeModel.DEFAULT_COVERAGE, 0, Integer.MAX_VALUE);
    	final SettingsModelBoolean ifmaxcoverage = new SettingsModelBoolean(Convert2AnnovarNodeModel.CFGKEY_IFMAXCOVERAGE, false);
    	final SettingsModelIntegerBounded maxcoverage = new SettingsModelIntegerBounded(Convert2AnnovarNodeModel.CFGKEY_MAXCOVERAGE, Convert2AnnovarNodeModel.DEFAULT_MAXCOVERAGE, 0, Integer.MAX_VALUE);
    	final SettingsModelBoolean includeinfo = new SettingsModelBoolean(Convert2AnnovarNodeModel.CFGKEY_INCLUDEINFO, false);
    	final SettingsModelBoolean allelicfrac = new SettingsModelBoolean(Convert2AnnovarNodeModel.CFGKEY_ALLELICFRAC, false);
    	final SettingsModelBoolean species = new SettingsModelBoolean(Convert2AnnovarNodeModel.CFGKEY_SPECIES, true);
    	final SettingsModelBoolean allallele = new SettingsModelBoolean(Convert2AnnovarNodeModel.CFGKEY_ALLALLELE, false);
    	final SettingsModelBoolean withzyg = new SettingsModelBoolean(Convert2AnnovarNodeModel.CFGKEY_WITHZYG, false);

    	createNewGroup("Annovar installpath.");
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(Convert2AnnovarNodeModel.CFGKEY_INSTALLPATH,null), "his_conv_id0", 0, true));
    	createNewGroup("Input file [possible formats: vcf, pileup, tsv or gff].");
    	addDialogComponent(new DialogComponentFileChooser(inFile, "his0_conv_id1", 0, "vcf","pileup","tsv","gff"));
    	createNewGroup("General options");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(ifmaxcoverage, "Use maximum coverage threshold:"));
    	addDialogComponent(new DialogComponentNumber(maxcoverage,"", /*step*/ 1));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentBoolean(includeinfo, "Include supporting information in output."));
    	createNewGroup("Options for pileup file format");
    	addDialogComponent(new DialogComponentNumber(snpqual,"Quality score threshold:", /*step*/ 1));
    	addDialogComponent(new DialogComponentNumber(coverage,"Read coverage threshold:", /*step*/ 1));
    	addDialogComponent(new DialogComponentBoolean(allelicfrac, "Print out allelic fraction rather than het/hom status."));
    	createNewGroup("Options for GFF3-SOLiD file format");
    	addDialogComponent(new DialogComponentNumber(snppvalue,"SNP P-value threshold:", /*step*/ 0.0001));
    	addDialogComponent(new DialogComponentBoolean(species, "If human, convert chr23/24/25 to X/Y/M."));
    	createNewGroup("Options for VCF4 file format");
    	addDialogComponent(new DialogComponentBoolean(allallele, "Print all alleles when multiple calls are present."));
    	addDialogComponent(new DialogComponentBoolean(withzyg, "Print zygosity when 'includeinfo' is used."));

    	snpqual.setEnabled(false);
    	coverage.setEnabled(false);
    	allelicfrac.setEnabled(false);
    	species.setEnabled(false);
    	snppvalue.setEnabled(false);
    	allallele.setEnabled(false);
    	withzyg.setEnabled(false);
    	
    	inFile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				String path2inFile = inFile.getStringValue();
				String fileExtension = path2inFile.substring(path2inFile.lastIndexOf(".")+1,path2inFile.length());
	    		snpqual.setEnabled(false);
	        	coverage.setEnabled(false);
	        	allelicfrac.setEnabled(false);
	        	species.setEnabled(false);
	        	snppvalue.setEnabled(false);
	        	allallele.setEnabled(false);
	        	withzyg.setEnabled(false);
		    	if(fileExtension.equals("vcf")) {
		    		if(includeinfo.getBooleanValue()) {
		    			withzyg.setEnabled(true);
		    		}
		        	allallele.setEnabled(true);
		    	} else if(fileExtension.equals("pileup")) {
		    	//	snpqual.setEnabled(true);
		        	coverage.setEnabled(true);
		        	allelicfrac.setEnabled(true);
		    	} else if(fileExtension.equals("gff")) {
		    		species.setEnabled(true);
		        	snppvalue.setEnabled(true);
		    	}
			}
		});
    	
    	ifmaxcoverage.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				maxcoverage.setEnabled(ifmaxcoverage.getBooleanValue());
			}
		});
    	
    	includeinfo.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
			String path2inFile = inFile.getStringValue();
			String fileExtension = path2inFile.substring(path2inFile.lastIndexOf(".")+1,path2inFile.length());
				withzyg.setEnabled(includeinfo.getBooleanValue() && fileExtension.equals("vcf"));
			}
		});
    	
    	
    	
    }
}

