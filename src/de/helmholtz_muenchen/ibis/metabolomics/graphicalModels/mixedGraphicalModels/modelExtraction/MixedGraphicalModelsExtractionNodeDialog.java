package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.modelExtraction;

import org.knime.core.data.DataTableSpec;
import org.knime.core.node.NodeSettingsRO;
import org.knime.core.node.NotConfigurableException;
import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;

/**
 * 
 * @author Jonas Zierer
 */
public class MixedGraphicalModelsExtractionNodeDialog extends DefaultNodeSettingsPane {

	/** SETTING MODELS */
	private final SettingsModelDoubleBounded  m_percIncl = new SettingsModelDoubleBounded(MixedGraphicalModelsExtractionNodeModel.CFGKEY_PERC_INCLUSION, 0.8, 0.0, 1.0);
	// FWER control only
	private final SettingsModelIntegerBounded m_ev       = new SettingsModelIntegerBounded(MixedGraphicalModelsExtractionNodeModel.CFGKEY_EV, 5, 0, Integer.MAX_VALUE);
	// Empirical p-values only
	private final SettingsModelIntegerBounded m_nEdges   = new SettingsModelIntegerBounded(MixedGraphicalModelsExtractionNodeModel.CFGKEY_NUM_OF_EDGES, 5, 1, Integer.MAX_VALUE);
	// Cores
	private final SettingsModelIntegerBounded m_cores    = new SettingsModelIntegerBounded(MixedGraphicalModelsExtractionNodeModel.CFGKEY_CORES, 3, 0, Integer.MAX_VALUE);

    
	protected MixedGraphicalModelsExtractionNodeDialog() {
        super();

        addDialogComponent(new DialogComponentNumber(
        		m_percIncl, 
        		"Inclusion Percentage", /* step size */0.1)
        );
        
        
        // FWER CONTROL ONLY
        this.createNewGroup("Intrinsic False-Positive Control (FWER)");
        addDialogComponent(new DialogComponentNumber(
        		m_ev,
        		"E(V):", /*step*/ 1, /*componentwidth*/ 5)
        );
        this.closeCurrentGroup();
        
        
        // EMPIRICAL P_VALUES ONLY
        this.createNewGroup("Number of edges per model");
        addDialogComponent(new DialogComponentNumber(
        		m_nEdges,
        		"Number of Edges selected from each model", /*step*/ 1, /*componentwidth*/ 5)
        );
        this.closeCurrentGroup();
        
        // PARALLEL EXECUTION
        this.createNewGroup("Parallelization");
        addDialogComponent(new DialogComponentNumber(
        		m_cores,
        		"Number of Cores", /*step*/ 1, /*componentwidth*/ 5)
        );
        this.closeCurrentGroup();
    }
	
	
	
	@Override
	public void loadAdditionalSettingsFrom(final NodeSettingsRO settings, final DataTableSpec[] specs) throws NotConfigurableException {
		// optional input exists --> Empirical P-Values	
		if(specs[1] != null & specs[1].getNumColumns()!=0){
			m_nEdges.setEnabled(true);
			m_ev.setEnabled(false);
		}else{
			m_nEdges.setEnabled(false);
			m_ev.setEnabled(true);
		}
		
	}
}

