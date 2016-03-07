package de.helmholtz_muenchen.ibis.ngs.Bcl2FastQ;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;

import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;


import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "Bcl2FastQ" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Kaarin Ahomaa
 */
public class Bcl2FastQNodeDialog extends HTExecutorNodeDialog {

	private final SettingsModelString ToolPath = new SettingsModelString(Bcl2FastQNodeModel.CFGKEY_TOOL_PATH,"");
	private final SettingsModelString InputPath = new SettingsModelString(Bcl2FastQNodeModel.CFGKEY_INPUT_PATH,"");
	private final SettingsModelString outfolder = new SettingsModelString(Bcl2FastQNodeModel.CFGKEY_OUTFOLDER,"");
	private final SettingsModelBoolean IsPaired = new SettingsModelBoolean(Bcl2FastQNodeModel.CFGKEY_ISPAIRED, true);
	private final SettingsModelIntegerBounded threads = new SettingsModelIntegerBounded(Bcl2FastQNodeModel.CFGKEY_THREADS,4, 1, Integer.MAX_VALUE);
	//private final SettingsModelString interop = new SettingsModelString(Bcl2FastQNodeModel.CFGKEY_INTEROP,"");
	
	
	
    /**
     * New pane for configuring the Bcl2FastQ node.
     */
    protected Bcl2FastQNodeDialog() {

    createNewGroup("Path to Tool");
    addDialogComponent(new DialogComponentFileChooser(ToolPath, "his_id_TOOL", 0));
    
    createNewGroup("Path to input files");
    addDialogComponent(new DialogComponentFileChooser(InputPath, "his_id_INPUT", 0, true));
    
    createNewGroup("Path to output files");
    addDialogComponent(new DialogComponentFileChooser(outfolder, "his_id_OUTPUT", 0, true));
    
    createNewGroup("Further options");
    addDialogComponent(new DialogComponentBoolean(IsPaired, "Paired end reads?"));
    addDialogComponent(new DialogComponentNumber(threads, "Number of threads", 1));
    
    //createNewGroup("Path to interop");
    //addDialogComponent(new DialogComponentFileChooser(interop, "his_id_INTEROP", 0, true));
    
    
    
   }
	@Override
	protected void updatePrefs() {
		// TODO Auto-generated method stub
		
	}
}

