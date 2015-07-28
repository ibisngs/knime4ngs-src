package de.helmholtz_muenchen.ibis.ngs.rawreadmanipulator;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;


/**
 * <code>NodeDialog</code> for the "RawReadManipulator" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class RawReadManipulatorNodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the RawReadManipulator node.
     */
 protected RawReadManipulatorNodeDialog() {
	
	 	final SettingsModelBoolean removeadapters = new SettingsModelBoolean(RawReadManipulatorNodeModel.CFGKEY_REMOVEADAPTERS, false);
		final SettingsModelString adapters = new SettingsModelString(RawReadManipulatorNodeModel.CFGKEY_ADAPTERS,"");
		final SettingsModelBoolean dotrimpolyat = new SettingsModelBoolean(RawReadManipulatorNodeModel.CFGKEY_DOTRIMPOLYAT, false);
    	final SettingsModelString trimpolyat = new SettingsModelString(RawReadManipulatorNodeModel.CFGKEY_TRIMPOLYAT,"");
	 	final SettingsModelBoolean filterfileexists = new SettingsModelBoolean(RawReadManipulatorNodeModel.CFGKEY_FILTERFILEEXISTS, true);
	 	final SettingsModelBoolean lengthcutoff = new SettingsModelBoolean(RawReadManipulatorNodeModel.CFGKEY_LENGTHCUTOFF, false);
    	final SettingsModelIntegerBounded minlength = new SettingsModelIntegerBounded(RawReadManipulatorNodeModel.CFGKEY_MINLENGTH, 30, 1, Integer.MAX_VALUE);
		final SettingsModelString preserve = new SettingsModelString(RawReadManipulatorNodeModel.CFGKEY_PRESERVE,"No");
		final SettingsModelBoolean ifillumina = new SettingsModelBoolean(RawReadManipulatorNodeModel.CFGKEY_IFILLUMINA, false);
		final SettingsModelString isillumina = new SettingsModelString(RawReadManipulatorNodeModel.CFGKEY_ISILLUMINA,"");
		final SettingsModelString convtophred = new SettingsModelString(RawReadManipulatorNodeModel.CFGKEY_CONVTOPHRED,"");
	 	final SettingsModelBoolean ifbarcodefile = new SettingsModelBoolean(RawReadManipulatorNodeModel.CFGKEY_IFBARCODEFILE, false);
		final SettingsModelString barcodefile = new SettingsModelString(RawReadManipulatorNodeModel.CFGKEY_BARCODEFILE,"");
		final SettingsModelBoolean usequalthreshold = new SettingsModelBoolean(RawReadManipulatorNodeModel.CFGKEY_USEQUALTHRESHOLD, false);
    	final SettingsModelIntegerBounded qualthreshold = new SettingsModelIntegerBounded(RawReadManipulatorNodeModel.CFGKEY_QUALTHRESHOLD, 1, 1, Integer.MAX_VALUE);
		final SettingsModelBoolean usetrimbyqual = new SettingsModelBoolean(RawReadManipulatorNodeModel.CFGKEY_USETRIMBYQUAL, false);
    	final SettingsModelIntegerBounded trimbyqual = new SettingsModelIntegerBounded(RawReadManipulatorNodeModel.CFGKEY_TRIMBYQUAL, 1, 1, Integer.MAX_VALUE);
	 	final SettingsModelBoolean useOtherFilterFile = new SettingsModelBoolean(RawReadManipulatorNodeModel.CFGKEY_USEOTHERFILTERFILE, false);
		final SettingsModelString otherFilterSettingsfile = new SettingsModelString(RawReadManipulatorNodeModel.CFGKEY_OTHERFILTERSETTINGSFILE,null);
		final SettingsModelBoolean trimbothends = new SettingsModelBoolean(RawReadManipulatorNodeModel.CFGKEY_TRIMBOTHENDS, false);
		final SettingsModelOptionalString outputfolder = new SettingsModelOptionalString(RawReadManipulatorNodeModel.CFGKEY_OUTPUTFOLDER, "", false);
		
    	/*adapters.setEnabled(false);
		trimpolyat.setEnabled(false);
		minlength.setEnabled(false);
		preserve.setEnabled(true);
		isillumina.setEnabled(false);
		convtophred.setEnabled(false);
		barcodefile.setEnabled(false);
		qualthreshold.setEnabled(false);
		trimbyqual.setEnabled(false);
		otherFilterSettingsfile.setEnabled(false);*/
    	
    	
    	createNewGroup("Filter-settings file (from modified FastQC)");
    	addDialogComponent(new DialogComponentBoolean(filterfileexists, "Use filter-settings calculated by FastQC (in-port)"));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(useOtherFilterFile, "Specify own filter-settings file"));
    	addDialogComponent(new DialogComponentFileChooser(otherFilterSettingsfile, "t1", 0, "filterSettings"));
    	
    	setHorizontalPlacement(false);
    	
    	createNewGroup("Barcode file (format: ID   BARCODE)");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(ifbarcodefile, "Split&trim by barcodes"));
    	addDialogComponent(new DialogComponentFileChooser(barcodefile, "Barcode file", 0, false));
    	setHorizontalPlacement(false);
    	
		createNewTab("Parameters");
    	
    	createNewGroup("Remove Adapters?");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(removeadapters, "Remove adapters"));
    	addDialogComponent(new DialogComponentString(adapters, "Adapter list (comma sep.)"));
    	setHorizontalPlacement(false);
    	
    	createNewGroup("Trim 3'/5' Poly-A/T prefixes:");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(dotrimpolyat, "Trim Poly-A/T"));
    	addDialogComponent(new DialogComponentStringSelection(trimpolyat, "Trim which end?","3\'","5\'"));
    	setHorizontalPlacement(false);
    	
    	createNewGroup("Discard reads with avg. quality below threshold");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usequalthreshold, "Set min. avg. quality"));
      	addDialogComponent(new DialogComponentNumber(
    			qualthreshold, "Min. avg. quality:", /*step*/ 1));
    	setHorizontalPlacement(false);
    	
    	createNewGroup("Trim reads from right until quality threshold is reached");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(usetrimbyqual, "Trim reads by quality"));
      	addDialogComponent(new DialogComponentNumber(
    			trimbyqual, "Quality threshold:", /*step*/ 1));
    	setHorizontalPlacement(false);
    	addDialogComponent(new DialogComponentBoolean(trimbothends, "Also trim from left"));	

    	
    	createNewGroup("Reads containing Ns:");
    	addDialogComponent(new DialogComponentStringSelection(
    			new SettingsModelString(RawReadManipulatorNodeModel.CFGKEY_REMOVEN,"No"), "Remove reads containing Ns?","No", "Yes"));
    	
    	createNewGroup("Keep single ends if partner is discarded?");
    	addDialogComponent(new DialogComponentStringSelection(preserve, "Keep single ends?","No", "Yes"));
    	
    	createNewGroup("Discard reads shorter than a given length (after processing):");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(lengthcutoff, "Set min. length"));
      	addDialogComponent(new DialogComponentNumber(
    			minlength, "Min. length", /*step*/ 1));
    	setHorizontalPlacement(false);
    	
    	createNewGroup("Options for Illumina-data:");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(ifillumina, "Data is from Illumina"));
      	addDialogComponent(new DialogComponentStringSelection(isillumina, "Choose quality encoding:","auto", "<1.3","1.3","1.5","1.9"));
    	setHorizontalPlacement(false);
      	addDialogComponent(new DialogComponentStringSelection(convtophred, "Convert quality to phred-based qualities?","No", "Yes"));
    	
    	createNewGroup("Thread number");
    	addDialogComponent(new DialogComponentNumber(
    			new SettingsModelIntegerBounded(
    				RawReadManipulatorNodeModel.CFGKEY_THREADCOUNT, 4, 1, Integer.MAX_VALUE),
    				"Number of threads to use:", /*step*/ 1));
    	
    	createNewGroup("Optional output folder for RRM");
    	addDialogComponent(new DialogComponentOptionalString(outputfolder, "output folder:"));
    	
    	
    	/**Checkboxes for the optional arguments**/
    	//checkbox for adapter removal
		 removeadapters.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent e) {
					adapters.setEnabled(removeadapters.getBooleanValue());
				}
		 });
		 //checkbox for polyAT trimming
		 dotrimpolyat.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent e) {
					trimpolyat.setEnabled(dotrimpolyat.getBooleanValue());
				}
		 });
		 //checkbox for Barcode File
		 
		 ifbarcodefile.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent e) {
					barcodefile.setEnabled(ifbarcodefile.getBooleanValue());
				}
		 });
		 //checkbox for length-cutoff
		 lengthcutoff.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent e) {
					minlength.setEnabled(lengthcutoff.getBooleanValue());
				}
		 });
		 //checkbox for Illumina quality encoding
		 ifillumina.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent e) {
					isillumina.setEnabled(ifillumina.getBooleanValue());
					convtophred.setEnabled(ifillumina.getBooleanValue());
				}
		 });
		 //checkbox for Trim by quality
		 usetrimbyqual.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent e) {
					trimbyqual.setEnabled(usetrimbyqual.getBooleanValue());
					trimbothends.setEnabled(usetrimbyqual.getBooleanValue());
					if(!usetrimbyqual.getBooleanValue()){
						trimbothends.setBooleanValue(false);
					}
				}
		 });
		 //checkbox for average read quality
		 usequalthreshold.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent e) {
					qualthreshold.setEnabled(usequalthreshold.getBooleanValue());
				}
		 });
		 
		 useOtherFilterFile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(useOtherFilterFile.getBooleanValue()){
					otherFilterSettingsfile.setEnabled(true);
					filterfileexists.setEnabled(false);
					filterfileexists.setBooleanValue(false);
					
				}
				else{
					otherFilterSettingsfile.setEnabled(false);
					useOtherFilterFile.setBooleanValue(false);
					filterfileexists.setEnabled(true);
				}
				

			}
		});
		 
		 filterfileexists.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(filterfileexists.getBooleanValue()){
					useOtherFilterFile.setEnabled(false);
					useOtherFilterFile.setBooleanValue(false);
					otherFilterSettingsfile.setEnabled(false);
				}
				else{
					useOtherFilterFile.setEnabled(true);
				}
				
			}
		});
		 

    }
}

