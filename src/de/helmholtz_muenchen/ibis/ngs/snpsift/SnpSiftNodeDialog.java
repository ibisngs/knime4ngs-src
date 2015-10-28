package de.helmholtz_muenchen.ibis.ngs.snpsift;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "SnpSift" Node.
 * 
 * @author Maximilian Hastreiter
 */
public class SnpSiftNodeDialog extends HTExecutorNodeDialog {

	private final SettingsModelString snpeff_folder = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_SNPEFF_FOLDER,"");
	private final SettingsModelString invcf = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_INVCF,"");
	private final SettingsModelString method = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_METHOD,"");
	
	/**Filter**/
	private final SettingsModelString filterstring = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_FILTERSTRING,"");
	private final SettingsModelDoubleBounded filterqual = new SettingsModelDoubleBounded(SnpSiftNodeModel.CFGKEY_FILTERQUAL, 20, 0, Double.MAX_VALUE);
	private final SettingsModelIntegerBounded filtercoverage = new SettingsModelIntegerBounded(SnpSiftNodeModel.CFGKEY_FILTERCOVERAGE, 10, 0, Integer.MAX_VALUE);
	private final SettingsModelBoolean filterqualbool = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_FILTERCOVERAGEBOOL, false);
	private final SettingsModelBoolean filtercoveragebool = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_FILTERCOVERAGEBOOL, false);
	
	/**Annotate**/
	private final SettingsModelOptionalString anninfo = new SettingsModelOptionalString(
			SnpSiftNodeModel.CFGKEY_ANNINFO,"",false);
	private final SettingsModelBoolean annid = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_ANNID, false);
	private final SettingsModelString annvcfdb = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_ANNVCFDB,"");
	
	/**TvTs**/
	private final SettingsModelString tstvhom = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_TSTVHOM,"");
	
	
	/**Intervals**/
	private final SettingsModelString interbed = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_INTERBED,"");
	private final SettingsModelBoolean interx = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_INTERX, false);

	/**dbnsfp**/
	private final SettingsModelString dbnsfp = new SettingsModelString(
			SnpSiftNodeModel.CFGKEY_DBNSFP,"");
	private final SettingsModelOptionalString dbnsfpfields = new SettingsModelOptionalString(
			SnpSiftNodeModel.CFGKEY_DBNSFPFFIELDS,"",false);
	private final SettingsModelBoolean dbnsfpfieldsall = new SettingsModelBoolean(SnpSiftNodeModel.CFGKEY_DBNSFPFFIELDSALL, false);

	
    /**
     * New pane for configuring the SnpSift node.
     */
    protected SnpSiftNodeDialog() {

    	super();
    	
    	createNewGroup("snpEff/snpSIFT directory");
    	addDialogComponent(new DialogComponentFileChooser(snpeff_folder, "par_1", 0, true));
    	createNewGroup("Input vcf file");
    	addDialogComponent(new DialogComponentFileChooser(invcf, "par_2", 0, false, ".vcf"));
   
    	addDialogComponent(new DialogComponentStringSelection(method, "Select method","Filter","Annotate","TsTv","Intervals","Annotate with dbnsfp"));
    	
    	createNewTab("Filter");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(filterqualbool, ""));
    	addDialogComponent(new DialogComponentNumber(filterqual, "Discard variants below QUAL value of:", 1));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(filtercoveragebool, ""));
    	addDialogComponent(new DialogComponentNumber(filtercoverage, "Discard variants below coverage value of:", 1));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentString(filterstring, "Enter additional filter criteria (SnpSift Syntax)"));
    	
    	filtercoverage.setEnabled(false);
    	filterqual.setEnabled(false);
    	setEnabled(false, "Filter");
		annvcfdb.setEnabled(false);
    	
    	
    	createNewTab("Annotate");
    	createNewGroup("Annotate using fields from a given vcf(e.g. dbSnp or 1000 genomes projects)");
    	addDialogComponent(new DialogComponentFileChooser(annvcfdb, "par_3", 0, false, ".vcf"));
    	closeCurrentGroup();
    	addDialogComponent(new DialogComponentBoolean(annid, "Only add ID fields (no INFO fields)"));
    	addDialogComponent(new DialogComponentOptionalString(anninfo, "Annotate using a list of info fields (list is a comma separated list of fields). Default: ALL."));
    	
    	createNewTab("TsTv");
    	addDialogComponent(new DialogComponentStringSelection(tstvhom, "Use only genotype fields of the selected type","any","hom","het"));
    	
    	createNewTab("Intervals");
    	createNewGroup("Specify interval file");
    	addDialogComponent(new DialogComponentFileChooser(interbed, "par_4", 0, false, ".bed",".txt"));
    	addDialogComponent(new DialogComponentBoolean(interx, "Exclude VCF entries in intervals"));    	
    	
    	createNewTab("Annotate with dbnsfp");
    	createNewGroup("Specify dbnsfp database");
    	addDialogComponent(new DialogComponentFileChooser(dbnsfp, "par_5", 0, false,".txt"));
    	addDialogComponent(new DialogComponentBoolean(dbnsfpfieldsall, "Annotate all fields, even if the database has an empty value"));   
    	addDialogComponent(new DialogComponentOptionalString(dbnsfpfields, "Specify which fields are used for annotation by a comma separated list of field names"));

    	
    	
    	filtercoveragebool.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(filtercoveragebool.getBooleanValue()){
						filtercoverage.setEnabled(true);
					}else{
						filtercoverage.setEnabled(false);
					}
			}
		});
    	filterqualbool.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(filterqualbool.getBooleanValue()){
						filterqual.setEnabled(true);
					}else{
						filterqual.setEnabled(false);
					}
			}
		});
    	dbnsfpfieldsall.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(dbnsfpfieldsall.getBooleanValue()){
						dbnsfpfields.setEnabled(false);
						dbnsfpfields.setIsActive(false);
					}else{
						dbnsfpfields.setEnabled(true);
					}
			}
		});
    	
    	method.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(method.getStringValue().equals("Filter")){
						setEnabled(true, "Filter");
					}else{
						setEnabled(false, "Filter");
					}
					if(method.getStringValue().equals("Annotate")){
						setEnabled(true, "Annotate");
						annvcfdb.setEnabled(true);						
					}else{
						setEnabled(false, "Annotate");
						annvcfdb.setEnabled(false);
					}
					if(method.getStringValue().equals("TsTv")){
						setEnabled(true, "TsTv");
					}else{
						setEnabled(false, "TsTv");
					}
					if(method.getStringValue().equals("Intervals")){
						setEnabled(true, "Intervals");
						interbed.setEnabled(true);
					}else{
						setEnabled(false, "Intervals");
						interbed.setEnabled(false);
					}
					if(method.getStringValue().equals("Annotate with dbnsfp")){
						setEnabled(true, "Annotate with dbnsfp");
						dbnsfp.setEnabled(true);
					}else{
						setEnabled(false, "Annotate with dbnsfp");
						dbnsfp.setEnabled(false);
					}		
			}
		});
    }


	@Override
	protected void updatePrefs() {
		// TODO Auto-generated method stub
		
	}
}

