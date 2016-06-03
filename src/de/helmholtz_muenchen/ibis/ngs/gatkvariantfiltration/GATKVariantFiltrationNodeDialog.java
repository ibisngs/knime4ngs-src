package de.helmholtz_muenchen.ibis.ngs.gatkvariantfiltration;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentOptionalString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelOptionalString;


import de.helmholtz_muenchen.ibis.utils.abstractNodes.GATKNode.GATKNodeDialog;

/**
 * <code>NodeDialog</code> for the "GATKVariantFiltration" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Maximilian Hastreiter
 */
public class GATKVariantFiltrationNodeDialog extends GATKNodeDialog {

	@Override
	protected void addDialogComponent() {
//		SettingsModelOptionalString QUAL= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_QUAL, "<50.0",true);
		SettingsModelOptionalString QD= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_QD, "<2.0",false);
		SettingsModelOptionalString FS= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_FS, ">60.0",false);
		SettingsModelOptionalString MQ= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_MQ, "<40.0",false);
		SettingsModelOptionalString HS= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_HS, ">13.0",false);
		SettingsModelOptionalString MQR= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_MQR, "<-12.5",false);
		SettingsModelOptionalString RPR= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_RPR, "<-8.0",false);
		SettingsModelOptionalString INFOFilterString = new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_INFOFilterString, "",false);
		SettingsModelOptionalString INFOFilterName = new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_INFOFilterName, "GATKVariantFiltration",false);
		
		SettingsModelBoolean NOCALL = new SettingsModelBoolean(GATKVariantFiltrationNodeModel.CFGKEY_NOCALL, true);
		SettingsModelOptionalString DP= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_DP, "<8.0",true);
		SettingsModelOptionalString GQ= new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_GQ, "<20.0",true);
		SettingsModelOptionalString FORMATFilterString = new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_FORMATFilterString, "",false);
		SettingsModelOptionalString FORMATFilterName = new SettingsModelOptionalString(GATKVariantFiltrationNodeModel.CFGKEY_FORMATFilterName, "GATKVariantFiltration",false);
		
		setDefaultTabTitle("VariantFiltration");	
//createNewTab("VariantFiltration");
    	
    	createNewGroup("Filter INFO field (variants)");
    	
    	addDialogComponent(new DialogComponentOptionalString(INFOFilterName, "Filter name"));

//    	addDialogComponent(new DialogComponentOptionalString(QUAL, "Quality cutoff (QUAL)"));
    	
    	addDialogComponent(new DialogComponentOptionalString(QD, "Variant Confidence/Quality by Depth (QD)"));
    	addDialogComponent(new DialogComponentOptionalString(FS, "Strand Bias (FS)"));
    	addDialogComponent(new DialogComponentOptionalString(MQ, "RMS Mapping Quality (MQ)"));
    	addDialogComponent(new DialogComponentOptionalString(HS, "HaplotypeScore"));
    	addDialogComponent(new DialogComponentOptionalString(MQR, "MappingQualityRankSum"));
    	addDialogComponent(new DialogComponentOptionalString(RPR, "ReadPosRankSum"));
    	addDialogComponent(new DialogComponentOptionalString(INFOFilterString, "Additional Filter Options",30));
	
    	
    	createNewGroup("Filter FORMAT field (genotypes)");
    	addDialogComponent(new DialogComponentOptionalString(FORMATFilterName, "Filter name"));
    	
    	addDialogComponent(new DialogComponentOptionalString(DP, "Coverage Cutoff (DP)"));
    	addDialogComponent(new DialogComponentOptionalString(GQ, "Genotype Quality (GQ)"));
    	addDialogComponent(new DialogComponentOptionalString(FORMATFilterString, "Additional Filter Options",30));
    	addDialogComponent(new DialogComponentBoolean(NOCALL,"Set filtered variants to no-call?"));

	}
	
	
}

