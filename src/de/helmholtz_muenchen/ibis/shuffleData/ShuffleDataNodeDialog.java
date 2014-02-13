package de.helmholtz_muenchen.ibis.shuffleData;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentColumnFilter;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelFilterString;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;

/**
 * <code>NodeDialog</code> for the "ShuffleData" Node.
 * Shuffle data within columns of the input data matrix.
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Jonas Zierer
 */
public class ShuffleDataNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring the ShuffleData node.
     */
	protected ShuffleDataNodeDialog() {
    	super();
    	

        // INT
        addDialogComponent(new DialogComponentNumber(
        		new SettingsModelIntegerBounded(ShuffleDataNodeModel.CFGKEY_RANDOM_SEED, 5, Integer.MIN_VALUE, Integer.MAX_VALUE),
        		"Random Seed:", /*step*/ 1, /*componentwidth*/ 5)
        );

        // LIST
        this.createNewGroup("Column selection");
        addDialogComponent(new DialogComponentColumnFilter(
        		new SettingsModelFilterString(ShuffleDataNodeModel.CFGKEY_COLUMNS),
        		/* # of input file from which to choose column names*/ 0, true)
        );
        this.closeCurrentGroup();
    }
}

