package de.helmholtz_muenchen.ibis.ngs.vat;

import javax.swing.JFileChooser;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "VAT" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */
public class VATNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring VAT node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
	
	final SettingsModelString VATFolder = new SettingsModelString(VATNodeModel.CFGKEY_VAT_FOLDER, VATNodeModel.DEF_VAT_FOLDER);
	final SettingsModelString intervals = new SettingsModelString(VATNodeModel.CFGKEY_INTERVALS, VATNodeModel.DEF_INTERVALS);
	final SettingsModelString transcripts = new SettingsModelString(VATNodeModel.CFGKEY_TRANSCRIPTS, VATNodeModel.DEF_TRANSCRIPTS);

	
    protected VATNodeDialog() {
        super();
        
        createNewGroup("Folder containing VAT executables");
        addDialogComponent(new DialogComponentFileChooser(VATFolder, "vatfolder", JFileChooser.OPEN_DIALOG, true));
        createNewGroup("Path to gene intervals");
        addDialogComponent(new DialogComponentFileChooser(intervals, "geneint", JFileChooser.OPEN_DIALOG, false, ".intervals|.interval"));
        createNewGroup("Path to transcript sequences");
        addDialogComponent(new DialogComponentFileChooser(transcripts, "transcripts", JFileChooser.OPEN_DIALOG, false, ".fa|.fasta"));

                    
    }
}

