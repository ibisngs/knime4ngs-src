package de.helmholtz_muenchen.ibis.ngs.thundercall;


import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentString;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
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

	private final SettingsModelString THUNDER = new SettingsModelString(ThunderCallNodeModel.CFGKEY_THUNDER_PATH, "");
	private final SettingsModelString BASE_NAME = new SettingsModelString(ThunderCallNodeModel.CFGKEY_BASE_NAME, "");
	private final SettingsModelDoubleBounded POST_PROB = new SettingsModelDoubleBounded(ThunderCallNodeModel.CFGKEY_POST_PROB, ThunderCallNodeModel.DEFAULT_POST_PROB, 0.1, 1.0);
	private final SettingsModelIntegerBounded MIN_DEPTH = new SettingsModelIntegerBounded(ThunderCallNodeModel.CFGKEY_MIN_DEPTH, ThunderCallNodeModel.DEFAULT_MIN_DEPTH, 1, 100);
    private final SettingsModelIntegerBounded MAX_DEPTH = new SettingsModelIntegerBounded(ThunderCallNodeModel.CFGKEY_MAX_DEPTH, ThunderCallNodeModel.DEFAULT_MAX_DEPTH, 100, 10000);

    protected ThunderCallNodeDialog() {
        super();
        
        createNewGroup("Path to thunder file");
    	DialogComponentFileChooser thunder= new DialogComponentFileChooser(THUNDER, "thunder", "");
    	addDialogComponent(thunder);
        @SuppressWarnings("deprecation")
		String thunderPath = BinaryHandler.checkToolAvailability("GPT_Freq");
    	if(thunderPath == null) {
    		thunderPath = "thunder GPT_Freq binary not found!";
    	}
    	THUNDER.setStringValue(thunderPath);
    	
    	addDialogComponent(new DialogComponentString(BASE_NAME, "Outfile Suffix (e.g. thunder10)"));

      	addDialogComponent(new DialogComponentNumber(POST_PROB,"Posterior probability:", /*step*/ 0.1, /*componentwidth*/ 5));

      	addDialogComponent(new DialogComponentNumber(MIN_DEPTH,"Minimum depth:", /*step*/ 10, /*componentwidth*/ 5));

      	addDialogComponent(new DialogComponentNumber(MAX_DEPTH,"Maximum depth:", /*step*/ 10, /*componentwidth*/ 5));

      	THUNDER.setEnabled(true);
        BASE_NAME.setEnabled(true);
        POST_PROB.setEnabled(true);
        MIN_DEPTH.setEnabled(true);
        MAX_DEPTH.setEnabled(true);
                    
    }
}

