package de.helmholtz_muenchen.ibis.ngs.bwa;

import java.io.File;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;
import org.knime.core.node.util.FilesHistoryPanel;

import de.helmholtz_muenchen.ibis.utils.BinaryHandler;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "BWA" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 */
public class BWANodeDialog extends HTExecutorNodeDialog {

	private final SettingsModelBoolean usePrefPage = new SettingsModelBoolean(BWANodeModel.CFGKEY_USEPREFPAGE,true);	
	private final SettingsModelString bwa = new SettingsModelString(BWANodeModel.CFGKEY_BWAFILE,"bwa");
	
	private final DialogComponentFileChooser dcfc = new DialogComponentFileChooser(bwa,"his_id_BWAPATH", 0);
	
	private final SettingsModelString readType = new SettingsModelString(BWANodeModel.CFGKEY_READTYPE,"auto-detect");
	private final SettingsModelString refseq = new SettingsModelString(BWANodeModel.CFGKEY_REFSEQFILE,null);
	private final SettingsModelBoolean indexrefseq = new SettingsModelBoolean(BWANodeModel.CFGKEY_CHECKINDEX, true);
	private final SettingsModelString readGroup = new SettingsModelString(BWANodeModel.CFGKEY_READGROUP, "@RG\\tID:foo\\tSM:bar");
	private final SettingsModelBoolean readGroupBoolean = new SettingsModelBoolean(BWANodeModel.CFGKEY_READGROUPBOOLEAN, false);
	private final SettingsModelIntegerBounded ALN_THREADS = new SettingsModelIntegerBounded(BWANodeModel.CFGKEY_THREADS,4, 0, Integer.MAX_VALUE);

	
    /**
     * New pane for configuring the BWA node.
     */
    protected BWANodeDialog() {
    	super();
    	    	
    	readGroup.setEnabled(false);
    	
    	createNewGroup("BWA");
    	
    	addDialogComponent(new DialogComponentBoolean(usePrefPage,"Use vales from KNIME4NGS preference page?"));
    	
    	addDialogComponent(dcfc);
    	dcfc.setEnabled(false);
    	createNewGroup("Reference sequence: FastA file (e.g. genome)");
    	addDialogComponent(new DialogComponentFileChooser(refseq, "his1_id_BWA", 0, ""));
    	createNewGroup("Options");
    	addDialogComponent(new DialogComponentBoolean(indexrefseq, "Index reference sequence (Has to be done if index does not exist yet)."));
    	addDialogComponent(new DialogComponentStringSelection(new SettingsModelString(BWANodeModel.CFGKEY_BWTINDEX,"BWT-SW"),"Algorithm for constructing BWT index:","BWT-SW","IS"));
    	addDialogComponent(new DialogComponentStringSelection(new SettingsModelString(BWANodeModel.CFGKEY_ALNALGO,"BWA-backtrack"),"Algorithm for mapping:","BWA-backtrack","BWA-SW","BWA-MEM"));
    	addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean(BWANodeModel.CFGKEY_CHECKCOLORSPACED, false), "Build color-space index."));
    	addDialogComponent(new DialogComponentStringSelection(readType,"Type of reads/ mapping:","auto-detect","single-end","paired-end"));
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(readGroupBoolean,"Specify Read Group Header:"));
    	addDialogComponent(new DialogComponentString(readGroup,""));
    	setHorizontalPlacement(false);
    	
    	addDialogComponent(new DialogComponentNumber(ALN_THREADS, "Number of Threads for BWA aln", 1));
    	
    	readType.setEnabled(false);
    	refseq.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
    	    	String path2refSeq = refseq.getStringValue();
    			File file0 = new File(path2refSeq+".amb");
    			File file1 = new File(path2refSeq+".ann");
    			File file2 = new File(path2refSeq+".bwt");
    			File file3 = new File(path2refSeq+".pac");
    			File file4 = new File(path2refSeq+".sa");
    			if(file0.exists() && file1.exists() && file2.exists() && file3.exists() && file4.exists()){
    				indexrefseq.setBooleanValue(false);
    				indexrefseq.setEnabled(true);
    			} else {
    				indexrefseq.setBooleanValue(true);
    				indexrefseq.setEnabled(false);
    			}
			}
		});
    	
    	readGroupBoolean.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					if(readGroupBoolean.getBooleanValue()){
						readGroup.setEnabled(true);
					}else{
						readGroup.setEnabled(false);
					}
			}
		});
    	
    	usePrefPage.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent arg0) {
				dcfc.setEnabled(!usePrefPage.getBooleanValue());
				if(usePrefPage.getBooleanValue()) {
					String toolPath = BinaryHandler.checkToolAvailability("bwa");
			    	if(toolPath == null) {
			    		toolPath = "";
			    	}
			    	bwa.setStringValue(toolPath);
				}
			}
    	});
    }
    
   
    public void onOpen() {
    	if(usePrefPage.getBooleanValue()){
	    	String toolPath = BinaryHandler.checkToolAvailability("bwa");
	    	if(toolPath == null) {
	    		toolPath = "";
	    	}
	    	bwa.setStringValue(toolPath);
    	}
    }
    
			
}

