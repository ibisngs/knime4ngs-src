package de.helmholtz_muenchen.ibis.ngs.vcftoolsfilter;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelInteger;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.knime.IBISKNIMENodesPlugin;
import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;


/**
 * <code>NodeDialog</code> for the "VEPFilter" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author tim.jeske
 */
public class VCFToolsFilterNodeDialog extends HTExecutorNodeDialog {

    /**
     * New pane for configuring the VCFFilter node.
     */

	private final SettingsModelString vcf_tools = new SettingsModelString(VCFToolsFilterNodeModel.CFGKEY_VCFTOOLS,"-");
	private final SettingsModelBoolean fill_an_ac = new SettingsModelBoolean(VCFToolsFilterNodeModel.CFGKEY_FILL_AN_AC,false);

	
	//genotype filter
	private final SettingsModelBoolean filter_by_DP = new SettingsModelBoolean(VCFToolsFilterNodeModel.CFGKEY_FILTER_BY_DP,false);
	private final SettingsModelInteger DP_threshold = new SettingsModelInteger(VCFToolsFilterNodeModel.CFGKEY_DP_THRESHOLD,8);
	
	private final SettingsModelBoolean filter_by_GQ = new SettingsModelBoolean(VCFToolsFilterNodeModel.CFGKEY_FILTER_BY_GQ,false);
	private final SettingsModelIntegerBounded GQ_threshold = new SettingsModelIntegerBounded(VCFToolsFilterNodeModel.CFGKEY_GQ_THRESHOLD,20,0,99);
	
	private final SettingsModelBoolean filter_by_AD = new SettingsModelBoolean(VCFToolsFilterNodeModel.CFGKEY_FILTER_BY_AD,false);
	private final SettingsModelInteger AD_threshold = new SettingsModelInteger(VCFToolsFilterNodeModel.CFGKEY_AD_THRESHOLD,8);
	
	//variant filter
	private final SettingsModelBoolean filter_by_GQ_mean = new SettingsModelBoolean(VCFToolsFilterNodeModel.CFGKEY_FILTER_BY_GQ_MEAN,false);
	private final SettingsModelIntegerBounded GQ_MEAN_threshold = new SettingsModelIntegerBounded(VCFToolsFilterNodeModel.CFGKEY_GQ_MEAN_THRESHOLD,35,0,99);
	
	private final SettingsModelBoolean filter_pass = new SettingsModelBoolean(VCFToolsFilterNodeModel.CFGKEY_FILTER_PASS,false);
	
	private final SettingsModelBoolean filter_by_callRate = new SettingsModelBoolean(VCFToolsFilterNodeModel.CFGKEY_FILTER_CALL_RATE,false);
	private final SettingsModelDoubleBounded callRate_threshold = new SettingsModelDoubleBounded(VCFToolsFilterNodeModel.CFGKEY_CALL_RATE_THRESHOLD,0.88,0.0,1.0);
		
    protected VCFToolsFilterNodeDialog() {
    	
    	addPrefPageSetting(vcf_tools, IBISKNIMENodesPlugin.VCFTOOLS);
    	
    	createNewGroup("Path to VCFtools binary");
    	addDialogComponent(new DialogComponentFileChooser(vcf_tools, "his_id_vcftools", 0, ""));
    	
    	createNewGroup("Further options");
    	addDialogComponent(new DialogComponentBoolean(fill_an_ac,"Restore allele counts?"));

    	
    	//genotype filter
    	createNewTab("Genotype Filter");
		addDialogComponent(new DialogComponentBoolean(filter_by_DP, "Filter genotypes by DP?"));
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentNumber(DP_threshold, "DP threshold",1));
		setHorizontalPlacement(false);
		
		addDialogComponent(new DialogComponentBoolean(filter_by_AD, "Filter genotypes by AD? (Pindel VCF)"));
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentNumber(AD_threshold, "AD threshold",1));
		setHorizontalPlacement(false);
		
		addDialogComponent(new DialogComponentBoolean(filter_by_GQ, "Filter genotypes by GQ?"));
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentNumber(GQ_threshold, "GQ threshold",1));
		setHorizontalPlacement(false);
		
		//variant filter
		createNewTab("Variant Filter");
		addDialogComponent(new DialogComponentBoolean(filter_pass,"Filter by FILTER is PASS?"));
	
		addDialogComponent(new DialogComponentBoolean(filter_by_callRate,"Filter by call rate?"));
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentNumber(callRate_threshold,"Call rate threshold",0.01));
		setHorizontalPlacement(false);
		
		addDialogComponent(new DialogComponentBoolean(filter_by_GQ_mean, "Filter genotypes by mean GQ?"));
		setHorizontalPlacement(true);
		addDialogComponent(new DialogComponentNumber(GQ_MEAN_threshold, "Mean GQ threshold",1));
		setHorizontalPlacement(false);
    }

    
//	@Override
//	protected void updatePrefs() {
//		if(usePrefPage.getBooleanValue()) {
//	    	String vcftoolsPath = IBISKNIMENodesPlugin.getDefault().getToolPathPreference("vcftools");
//	    	if(vcftoolsPath != null && !vcftoolsPath.equals("")) {
//	    		vcf_tools.setStringValue(vcftoolsPath);
//	    		vcf_tools.setEnabled(false);
//	    	} else {
//	    		vcf_tools.setEnabled(true);
//	    	}
//	    	
//		} else {
//			vcf_tools.setEnabled(true);
//		}
//	}
}