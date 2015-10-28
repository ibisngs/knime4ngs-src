package de.helmholtz_muenchen.ibis.ngs.snpeff;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTExecutorNodeDialog;

/**
 * <code>NodeDialog</code> for the "SnpEff" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author
 */
public class SnpEffNodeDialog extends HTExecutorNodeDialog {

	/*Mandatory options*/
	final SettingsModelString snpeff_folder = new SettingsModelString(
			SnpEffNodeModel.CFGKEY_SNPEFF_FOLDER, null);
	final SettingsModelString database = new SettingsModelString(
			SnpEffNodeModel.CFGKEY_DATABASE, null);
	final SettingsModelString vcf_file = new SettingsModelString(
			
			SnpEffNodeModel.CFGKEY_VCF_FILE, null);
	/*Sequence change filter options*/
	final SettingsModelBoolean useminq = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_USEMINQ, false);
	final SettingsModelDoubleBounded minq = new SettingsModelDoubleBounded(
			SnpEffNodeModel.CFGKEY_MINQ, 0.0, 0.0, Double.MAX_VALUE);
	final SettingsModelBoolean useminc = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_USEMINC, false);
	final SettingsModelIntegerBounded minc = new SettingsModelIntegerBounded(
			SnpEffNodeModel.CFGKEY_MINC, 1, 1, Integer.MAX_VALUE);
	
	final SettingsModelBoolean del = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_DEL, false);
	final SettingsModelBoolean ins = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_INS, false);
	final SettingsModelBoolean hom = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_HOM, false);
	final SettingsModelBoolean het = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_HET, false);
	final SettingsModelBoolean mnp = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_MNP, false);
	final SettingsModelBoolean snp = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_SNP, false);
	
	
	/*Results filter options*/
	final SettingsModelBoolean usebedfile = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_USEBEDFILE, false);
	final SettingsModelString bed_file = new SettingsModelString(
			SnpEffNodeModel.CFGKEY_BED_FILE, null);
	final SettingsModelBoolean no_downstream = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_DOWNSTREAM, false);
	final SettingsModelBoolean no_intergenic = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_INTERGENIC, false);
	final SettingsModelBoolean no_intronic = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_INTRONIC, false);
	final SettingsModelBoolean no_upstream = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_UPSTREAM, false);
	final SettingsModelBoolean no_utr = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_NO_UTR, false);
	
	/*Annotations options*/
	final SettingsModelBoolean lof = new SettingsModelBoolean(
			SnpEffNodeModel.CFGKEY_LOF, false);
	
    /**
     * New pane for configuring the SnpEff node.
     */
    protected SnpEffNodeDialog() {
    	
    	super();
    	
    	createNewGroup("snpEff directory");
    	addDialogComponent(new DialogComponentFileChooser(snpeff_folder, "par_1", 0, true));
    	
    	createNewGroup("Database name");
    	addDialogComponent(new DialogComponentString(database, ""));
    	
    	createNewGroup("Input vcf file");
    	addDialogComponent(new DialogComponentFileChooser(vcf_file, "par_2", 0, false, ".vcf"));
    	
    	createNewGroup("Annotation options");
    	addDialogComponent(new DialogComponentBoolean(lof, "Add loss of function (LOF) and Nonsense mediated decay (NMD) tags"));

    	
    	createNewTab("Sequence change filter");
    	createNewGroup("Set thresholds");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(useminq, "Set minimum quality value"));
    	addDialogComponent(new DialogComponentNumber(minq, "Threshold", /*step*/ 0.1));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(useminc, "Set minimum coverage"));
    	addDialogComponent(new DialogComponentNumber(minc,	"Threshold", /*step*/ 1));
    	setHorizontalPlacement(false);
    	createNewGroup("Filter by mutation type");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(del, "Analyze deletions only"));
    	addDialogComponent(new DialogComponentBoolean(ins, "Analyze insertions only"));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(hom, "Analyze homozygous variants only"));
    	addDialogComponent(new DialogComponentBoolean(het, "Analyze heterozygous variants only"));
    	setHorizontalPlacement(false);
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(snp, "Analyze SNPs only"));
    	addDialogComponent(new DialogComponentBoolean(mnp, "Analyze MNPs only"));
    	setHorizontalPlacement(false);
    	
    	
    	createNewTab("Results filter");
    	createNewGroup("Restrict variants by chromosomal location");
    	addDialogComponent(new DialogComponentBoolean(usebedfile, "Only analyze changes that intersect specified intervals from bed file"));
    	addDialogComponent(new DialogComponentFileChooser(bed_file, "par_3", 0, false, ".bed"));
    	createNewGroup("Filter features");
    	addDialogComponent(new DialogComponentBoolean(no_downstream, "Don't show DOWNSTREAM changes"));
    	addDialogComponent(new DialogComponentBoolean(no_intergenic, "Don't show INTERGENIC changes"));
    	addDialogComponent(new DialogComponentBoolean(no_intronic, "Don't show INTRONIC changes"));
    	addDialogComponent(new DialogComponentBoolean(no_upstream, "Don't show UPSTREAM changes"));
    	addDialogComponent(new DialogComponentBoolean(no_utr, "Don't show UTR changes"));
    	
    	/*Sequence change filter options change listener*/
    	useminq.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(useminq.getBooleanValue()){
					minq.setEnabled(true);
				}
				else{
					minq.setEnabled(false);
				}
			}
		});
    	useminc.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(useminc.getBooleanValue()){
					minc.setEnabled(true);
				}
				else{
					minc.setEnabled(false);
				}
			}
		});
    	
    	del.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(del.getBooleanValue()){
					ins.setEnabled(false);
				}
				else{
					ins.setEnabled(true);
				}
			}
		});
    	ins.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(ins.getBooleanValue()){
					del.setEnabled(false);
				}
				else{
					del.setEnabled(true);
				}
			}
		});
    	hom.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(hom.getBooleanValue()){
					het.setEnabled(false);
				}
				else{
					het.setEnabled(true);
				}
			}
		});
    	het.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(het.getBooleanValue()){
					hom.setEnabled(false);
				}
				else{
					hom.setEnabled(true);
				}
			}
		});
    	mnp.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(mnp.getBooleanValue()){
					snp.setEnabled(false);
				}
				else{
					snp.setEnabled(true);
				}
			}
		});
    	snp.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(snp.getBooleanValue()){
					mnp.setEnabled(false);
				}
				else{
					mnp.setEnabled(true);
				}
			}
		});
    	
    	
    	/*Results filter options change listener*/
    	usebedfile.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				if(usebedfile.getBooleanValue()){
					bed_file.setEnabled(true);
				}
				else{
					bed_file.setEnabled(false);
				}
			}
		});
    	
    }

	@Override
	protected void updatePrefs() {
		// TODO Auto-generated method stub
		
	}
}

