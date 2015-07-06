package de.helmholtz_muenchen.ibis.ngs.thundercall;

import javax.swing.JFileChooser;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.BinaryHandler;

/**
 * <code>NodeDialog</code> for the "ThunderCall" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tanzeem Haque
 */
public class ThunderCallNodeDialog extends DefaultNodeSettingsPane {

	private final SettingsModelString SAMTOOLS_HYBRID = new SettingsModelString(ThunderCallNodeModel.CFGKEY_SAMTOOLS_HYBRID_PATH, "");
	private final SettingsModelString THUNDER = new SettingsModelString(ThunderCallNodeModel.CFGKEY_THUNDER_PATH, "");
	private final SettingsModelString REF_GENOME = new SettingsModelString(ThunderCallNodeModel.CFGKEY_REF_GENOME, "");
	private final SettingsModelString BASE_NAME = new SettingsModelString(ThunderCallNodeModel.CFGKEY_BASE_NAME, "");
	private final SettingsModelDoubleBounded POST_PROB = new SettingsModelDoubleBounded(ThunderCallNodeModel.CFGKEY_POST_PROB, ThunderCallNodeModel.DEFAULT_POST_PROB, 0.1, 1.0);

    protected ThunderCallNodeDialog() {
        super();
        
        createNewGroup("Path to samtools-hybrid file");
        String samHybridPath = BinaryHandler.checkToolAvailability("samtools-hybrid");
    	if(samHybridPath == null) {
    		samHybridPath = "samtools-hybrid binary not found!";
    	}
    	SAMTOOLS_HYBRID.setStringValue(samHybridPath);
    	
    	createNewGroup("Path to thunder file");
        String thunderPath = BinaryHandler.checkToolAvailability("GPT_Freq");
    	if(thunderPath == null) {
    		thunderPath = "thunder GPT_Freq binary not found!";
    	}
    	THUNDER.setStringValue(thunderPath);
    	
    	createNewGroup("Reference Genome");
    	DialogComponentFileChooser ref_genome= new DialogComponentFileChooser(REF_GENOME, "ref_genome_variant_filter", JFileChooser.OPEN_DIALOG, false, ".txt|.fa|.fasta");
    	ref_genome.setBorderTitle("Choose the reference genome");
    	addDialogComponent(ref_genome);
    	
    	addDialogComponent(new DialogComponentString(BASE_NAME, "Outfile Suffix"));

      	addDialogComponent(new DialogComponentNumber(POST_PROB,"Posterior probability:", /*step*/ 0.1, /*componentwidth*/ 5));

                    
    }
}

