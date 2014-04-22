package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.mixedGraphicalModels.modelExtraction;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;

/**
 * 
 * @author Jonas Zierer
 */
public class MixedGraphicalModelsExtractionNodeDialog extends DefaultNodeSettingsPane {


	protected MixedGraphicalModelsExtractionNodeDialog() {
        super();
        

        addDialogComponent(new DialogComponentNumber(
        		new SettingsModelIntegerBounded(MixedGraphicalModelsExtractionNodeModel.CFGKEY_EV, 5, Integer.MIN_VALUE, Integer.MAX_VALUE),
        		"E(V):", /*step*/ 1, /*componentwidth*/ 5)
        );
        addDialogComponent(new DialogComponentNumber(
        		new SettingsModelDoubleBounded(MixedGraphicalModelsExtractionNodeModel.CFGKEY_PERC_INCLUSION, 0.8, 0.0, 1.0), 
        		"Inclusion Percentage", /* step size */0.1)
        );


    }
}

