package de.helmholtz_muenchen.ibis.ngs.filtersummary;

import javax.swing.JFileChooser;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "FilterSummary" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author 
 */

/*
 * options
 * filter settings:
 * -min coverage
 * -min fraction of supporting reads
 * 
 * annotation settings:
 * -omim ids + path to gene_info file
 * -pathways + path to wikipathway file
 * -number of interactions + reactome file
 */

	
public class FilterSummaryNodeDialog extends DefaultNodeSettingsPane {
	
	final SettingsModelIntegerBounded min_coverage = new SettingsModelIntegerBounded(FilterSummaryNodeModel.CFGKEY_MIN_COVERAGE, FilterSummaryNodeModel.DEF_MIN_COVERAGE, FilterSummaryNodeModel.MIN_MIN_COVERAGE, FilterSummaryNodeModel.MAX_MIN_COVERAGE);
	final SettingsModelDoubleBounded supp_reads_frac = new SettingsModelDoubleBounded(FilterSummaryNodeModel.CFGKEY_SUPP_READS_FRAC, FilterSummaryNodeModel.DEF_SUPP_READS_FRAC, FilterSummaryNodeModel.MIN_SUPP_READS_FRAC, FilterSummaryNodeModel.MAX_SUPP_READS_FRAC);
	final SettingsModelBoolean write_stats = new SettingsModelBoolean(FilterSummaryNodeModel.CFGKEY_WIRTE_STATS, FilterSummaryNodeModel.DEF_WRITE_STATS);
	final SettingsModelBoolean ann_omim= new SettingsModelBoolean(FilterSummaryNodeModel.CFGKEY_ANN_OMIM, FilterSummaryNodeModel.DEF_ANN_OMIM);
	final SettingsModelString gene_info_file = new SettingsModelString(FilterSummaryNodeModel.CFGKEY_GENE_INFO_FILE, FilterSummaryNodeModel.DEF_GENE_INFO_FILE);
	final SettingsModelBoolean ann_pathways = new SettingsModelBoolean(FilterSummaryNodeModel.CFGKEY_ANN_PATHWAY, FilterSummaryNodeModel.DEF_ANN_PATHWAY);
	final SettingsModelString wikipathway_file = new SettingsModelString(FilterSummaryNodeModel.CFGKEY_WIKIPATHWAY_FILE, FilterSummaryNodeModel.DEF_WIKIPATHWAY_FILE);
	final SettingsModelBoolean ann_num_ppi = new SettingsModelBoolean(FilterSummaryNodeModel.CFGKEY_ANN_NUM_PPI, FilterSummaryNodeModel.DEF_ANN_NUM_PPI);
	final SettingsModelString reactome_file = new SettingsModelString(FilterSummaryNodeModel.CFGKEY_REACTOME_FILE, FilterSummaryNodeModel.DEF_REACTOME_FILE);
	final SettingsModelBoolean org_vcf_col= new SettingsModelBoolean(FilterSummaryNodeModel.CFGKEY_ORG_VCF_COL, FilterSummaryNodeModel.DEF_ORG_VCF_COL);
    /**
     * New pane for configuring FilterSummary node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
	
    protected FilterSummaryNodeDialog() {
        super();
        
        createNewGroup("Filtering");
        addDialogComponent(new DialogComponentNumber(min_coverage, "Minimum coverage", 1, 5));
        addDialogComponent(new DialogComponentNumber(supp_reads_frac, "Minimum fraction of supporting reads", 0.01, 5));
        
        createNewGroup("Output");
        addDialogComponent(new DialogComponentBoolean(write_stats, "Write summary file with variant counts"));
        addDialogComponent(new DialogComponentBoolean(org_vcf_col, "Add original vcf info and sample column to variant table"));
        
        createNewGroup("Annotation");
        
        addDialogComponent(new DialogComponentBoolean(ann_omim, "Annotate OMIM IDs from NCBI gene_info file"));
        DialogComponentFileChooser dc = new DialogComponentFileChooser(gene_info_file, "geneinfo_file", JFileChooser.OPEN_DIALOG, false, ".gene_info");
        dc.setBorderTitle("NCBI gene_info file");
        addDialogComponent(dc);
        gene_info_file.setEnabled(false);

        
        // activate file chooser
        ann_omim.addChangeListener(new ChangeListener() {
			
			@Override
			public void stateChanged(ChangeEvent e) {
				gene_info_file.setEnabled(ann_omim.getBooleanValue());
			}
		});
        
        addDialogComponent(new DialogComponentBoolean(ann_pathways, "Annotate pathways from Wikipathways"));
        DialogComponentFileChooser dc2 = new DialogComponentFileChooser(wikipathway_file, "wikipath_file", JFileChooser.OPEN_DIALOG, false);
        dc2.setBorderTitle("Wikipathway file");
        addDialogComponent(dc2);
        wikipathway_file.setEnabled(false);
        
        // activate file chooser
        ann_pathways.addChangeListener(new ChangeListener() {
			
			@Override
			public void stateChanged(ChangeEvent e) {
				wikipathway_file.setEnabled(ann_pathways.getBooleanValue());				
			}
		});
        
        addDialogComponent(new DialogComponentBoolean(ann_num_ppi, "Annotate number of interactions from Reactome"));
        DialogComponentFileChooser dc3 = new DialogComponentFileChooser(reactome_file, "react_file", JFileChooser.OPEN_DIALOG, false);
        dc3.setBorderTitle("Reactome interaction file");
        addDialogComponent(dc3);
        reactome_file.setEnabled(false);
        
        //activate file chooser
        ann_num_ppi.addChangeListener(new ChangeListener() {
			
			@Override
			public void stateChanged(ChangeEvent e) {
				reactome_file.setEnabled(ann_num_ppi.getBooleanValue());			
			}
		});
                    
    }
}

