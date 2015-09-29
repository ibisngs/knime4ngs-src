package de.helmholtz_muenchen.ibis.ngs.vcfnormalizer;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "VCFNormalizer" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author tim.jeske
 */
public class VCFNormalizerNodeDialog extends HTExecutorNodeDialog {

	private final SettingsModelString vt_path = new SettingsModelString(VCFNormalizerNodeModel.CFGKEY_VT,"");
	private final SettingsModelString vcf_in = new SettingsModelString(VCFNormalizerNodeModel.CFGKEY_VCF, "");
	private final SettingsModelString ref_genome = new SettingsModelString(VCFNormalizerNodeModel.CFGKEY_REF_GENOME,"");
	
    /**
     * New pane for configuring the VCFNormalizer node.
     */
    protected VCFNormalizerNodeDialog() {

    	createNewGroup("Path to vt directory");
    	addDialogComponent(new DialogComponentFileChooser(vt_path, "his_id_VT_PATH", 0, ""));
    	
    	createNewGroup("Path to input VCF file");
    	addDialogComponent(new DialogComponentFileChooser(vcf_in, "his_id_VCF_IN", 0, ".vcf|.vcf.gz"));
    	
    	createNewGroup("Path to reference genome");
    	addDialogComponent(new DialogComponentFileChooser(ref_genome, "his_id_REF_GENOME", 0, ".fa|.fasta"));
    }

	@Override
	protected void updatePrefs() {
		if(usePrefPage.getBooleanValue()) {
	    	String vtPath = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("vt");
	    	if(vtPath != null && !vtPath.equals("")) {
	    		vt_path.setStringValue(vtPath);
	    		vt_path.setEnabled(false);
			} else {
				vt_path.setEnabled(true);
			}
	    	String refGenome = IBISKNIMENodesPlugin.getDefault().getRefGenomePreference();
	    	if(refGenome != null && !refGenome.equals("")) {
	    		ref_genome.setStringValue(refGenome);
	    		ref_genome.setEnabled(false);
	    	} else {
	    		ref_genome.setEnabled(true);
	    	}
		} else {
			vt_path.setEnabled(true);
			ref_genome.setEnabled(true);
		}
		
	}
}