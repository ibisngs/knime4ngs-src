package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.edgeRankingFDRbackground;


import java.awt.GridLayout;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

import org.knime.core.data.DataTableSpec;
import org.knime.core.data.IntValue;
import org.knime.core.node.InvalidSettingsException;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NodeSettingsWO;
import org.knime.core.node.NotConfigurableException;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;

import de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.edgeRanking.MixedGraphicalModelsEdgeRankingNodeDialog;
import de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.edgeRanking.MixedGraphicalModelsEdgeRankingNodeModel;
import de.helmholtz_muenchen.ibis.utils.Global;


/**
 * <code>NodeDialog</code> for the "GraphicalModelsNode" Node.
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jonas Zierer
 */
public class MixedGraphicalModelsEdgeRankingBackgroundNodeDialog extends MixedGraphicalModelsEdgeRankingNodeDialog {

	private final SpinnerNumberModel m_sampleNum_FDR = new SpinnerNumberModel(MixedGraphicalModelsEdgeRankingBackgroundNodeModel.DEFAULT_SAMPLENUM_FDR, 1, 10000, 1);
	private final JSpinner stabSel_sampleNumFDR_spinner = new JSpinner(m_sampleNum_FDR);
	
	private final JComboBox<String> m_randomSeedsFDRCol;
	private String[] m_availableCols_fdrRSeed = new String[]{MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_AVAILABLE_COLS};

	
    /**
     * Creates a new dialog.
     */
    MixedGraphicalModelsEdgeRankingBackgroundNodeDialog() {
        super();
        
        // INITIALIZE PANEL FOR STABILITY SELECTION PARAMETERS
        final JPanel fdrPanel = new JPanel(new GridLayout(0,2, 5, 10));
        
        // num samples
        final JLabel stabSel_sampleNum_label     = new JLabel("Number of Samples");
        fdrPanel.add(stabSel_sampleNum_label);
        fdrPanel.add(stabSel_sampleNumFDR_spinner);
        
        
        // Random seeds col if optional input is given
        JLabel label_randomSeed = new JLabel("Random Seed");
        fdrPanel.add(label_randomSeed);
        m_randomSeedsFDRCol = new JComboBox<String>(m_availableCols_fdrRSeed);
        fdrPanel.add(m_randomSeedsFDRCol);
        
        super.addTab(" FDR Estimation ", fdrPanel);
    }
    

    
	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void loadSettingsFrom(final NodeSettingsRO settings, final DataTableSpec[] specs) throws NotConfigurableException {
		super.loadSettingsFrom(settings, specs);
		
		// additional settings
		m_sampleNum_FDR.setValue(settings.getInt(MixedGraphicalModelsEdgeRankingBackgroundNodeModel.CFGKEY_SS_SAMPLE_FDR, MixedGraphicalModelsEdgeRankingBackgroundNodeModel.DEFAULT_SAMPLENUM_FDR));
		
		
		if(specs[1] != null){
			DataTableSpec optionalInputSpecs = (DataTableSpec)specs[MixedGraphicalModelsEdgeRankingBackgroundNodeModel.INPORT_RANDOM_SEEDS_EMP_PVALS];
			this.m_availableCols_fdrRSeed = Global.getValidCols(optionalInputSpecs, IntValue.class);
		}else{
			m_availableCols_fdrRSeed = new String[]{MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_AVAILABLE_COLS};
		}
		DefaultComboBoxModel<String> model = new DefaultComboBoxModel<String>( m_availableCols_fdrRSeed );
		this.m_randomSeedsFDRCol.setModel(model);
		
		if(this.m_availableCols_fdrRSeed[0].equals(MixedGraphicalModelsEdgeRankingNodeModel.DEFAULT_AVAILABLE_COLS)){
			this.m_randomSeedsFDRCol.setEnabled(false);
			this.stabSel_sampleNumFDR_spinner.setEnabled(true);
		}else{
			this.m_randomSeedsFDRCol.setEnabled(true);
			this.stabSel_sampleNumFDR_spinner.setEnabled(false);
			this.m_randomSeedsFDRCol.setSelectedItem(settings.getString(MixedGraphicalModelsEdgeRankingBackgroundNodeModel.CFGKEY_RAND_SEED_EP_COL, ""));	
		}
	}

	/**
	 * {@inheritDoc}
	 */
	@Override
	protected void saveSettingsTo(final NodeSettingsWO settings) throws InvalidSettingsException {
		super.saveSettingsTo(settings);
		
		// additional settings
		settings.addInt(MixedGraphicalModelsEdgeRankingBackgroundNodeModel.CFGKEY_SS_SAMPLE_FDR, m_sampleNum_FDR.getNumber().intValue());
		settings.addString(MixedGraphicalModelsEdgeRankingBackgroundNodeModel.CFGKEY_RAND_SEED_EP_COL, (String)this.m_randomSeedsFDRCol.getSelectedItem());

	}
}

