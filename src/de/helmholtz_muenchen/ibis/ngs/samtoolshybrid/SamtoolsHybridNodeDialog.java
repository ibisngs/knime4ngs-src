package de.helmholtz_muenchen.ibis.ngs.samtoolshybrid;

import javax.swing.JFileChooser;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.ngs.thundercall.ThunderCallNodeModel;
import de.helmholtz_muenchen.ibis.utils.BinaryHandler;

/**
 * <code>NodeDialog</code> for the "SamtoolsHybrid" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tanzeem Haque
 */
public class SamtoolsHybridNodeDialog extends DefaultNodeSettingsPane {

	private final SettingsModelString SAMTOOLS_HYBRID = new SettingsModelString(SamtoolsHybridNodeModel.CFGKEY_SAMTOOLS_HYBRID_PATH, "");
	private final SettingsModelString REF_GENOME = new SettingsModelString(SamtoolsHybridNodeModel.CFGKEY_REF_GENOME, "");
	
	protected SamtoolsHybridNodeDialog() {
	        
		super();
	        
	        createNewGroup("Path to samtools-hybrid file");
	    	DialogComponentFileChooser sam= new DialogComponentFileChooser(SAMTOOLS_HYBRID, "samtools-hybrid", "");
	    	addDialogComponent(sam);

	        String samHybridPath = BinaryHandler.checkToolAvailability("samtools-hybrid");
	    	if(samHybridPath == null) {
	    		samHybridPath = "samtools-hybrid binary not found!";
	    	}
	    	SAMTOOLS_HYBRID.setStringValue(samHybridPath);
	    	
	    	createNewGroup("Reference Genome");
	    	DialogComponentFileChooser ref_genome= new DialogComponentFileChooser(REF_GENOME, "ref_genome", JFileChooser.OPEN_DIALOG, false, ".txt|.fa|.fasta");
	    	ref_genome.setBorderTitle("Choose the reference genome");
	    	addDialogComponent(ref_genome);
	    	
	                    
	    }
}

