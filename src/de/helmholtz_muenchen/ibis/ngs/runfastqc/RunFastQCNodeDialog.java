package de.helmholtz_muenchen.ibis.ngs.runfastqc;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "RunFastQC" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class RunFastQCNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the RunFastQC node.
     */
    protected RunFastQCNodeDialog() {
    	
    	final SettingsModelString reads1 = new SettingsModelString(RunFastQCNodeModel.CFGKEY_READSEQFILE,null);
    	final SettingsModelString reads2 = new SettingsModelString(RunFastQCNodeModel.CFGKEY_READSEQFILE2,null);
    	final SettingsModelString readType = new SettingsModelString(RunFastQCNodeModel.CFGKEY_READTYPE,"single-end");
    	reads2.setEnabled(false);

    	createNewGroup("Read sequences (FastQ, FastA, or BAM file)");
    	addDialogComponent(new DialogComponentFileChooser(reads1, "his0_runfastqc_id", 0, ""));
    	createNewGroup("Second file (FastQ or FastA) is optional (paired-end mapping)");
    	addDialogComponent(new DialogComponentFileChooser(reads2, "his01_runfastqc_id", 0, ""));
    	createNewGroup("Parameter");
    	addDialogComponent(new DialogComponentStringSelection(readType,"Type of reads/ mapping:","single-end","paired-end"));
    			
        	readType.addChangeListener(new ChangeListener(){
        		public void stateChanged(final ChangeEvent e) {
        	    	if(reads1.getStringValue().length()>0){
        	    		String f1 = reads1.getStringValue();
        	    		reads2.setEnabled(readType.getStringValue().equals("paired-end") && !f1.substring(f1.length()-3,f1.length()).equals("bam"));
        	    		}
        	    	}	
        	});
    		
        	reads1.addChangeListener(new ChangeListener() {
    			public void stateChanged(final ChangeEvent e) {
        	    	if(reads1.getStringValue().length()>0){
    				String f1 = reads1.getStringValue();
    				reads2.setEnabled(!f1.substring(f1.length()-3,f1.length()).equals("bam") && readType.getStringValue().equals("paired-end"));
        	    	}
        	    	}
    		});


    	
    }
}

