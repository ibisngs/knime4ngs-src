package de.helmholtz_muenchen.ibis.ngs.art;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "Art" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Syeda Tanzeem Haque
 */
public class ArtNodeDialog extends DefaultNodeSettingsPane {

    /**
     * New pane for configuring Art node dialog.
     * This is just a suggestion to demonstrate possible default dialog
     * components.
     */
	final SettingsModelIntegerBounded length = new SettingsModelIntegerBounded(ArtNodeModel.CFGKEY_LENGTH, ArtNodeModel.DEFAULT_LENGTH,1, 250);
    final SettingsModelIntegerBounded mean_size = new SettingsModelIntegerBounded(ArtNodeModel.CFGKEY_MEAN_SIZE, ArtNodeModel.DEFAULT_MEAN_SIZE,1, Integer.MAX_VALUE);
    final SettingsModelIntegerBounded fold = new SettingsModelIntegerBounded(ArtNodeModel.CFGKEY_FOLD, ArtNodeModel.DEFAULT_FOLD,1, 100);
    final SettingsModelIntegerBounded sd = new SettingsModelIntegerBounded(ArtNodeModel.CFGKEY_SD, ArtNodeModel.DEFAULT_LENGTH,1, 100);
//
//    final SettingsModelBoolean NO_MASK = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_NO_MASK, false); //if true, use -nf 0
//	final SettingsModelBoolean MATE_PAIR = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_MATE_PAIR, false); //if true, use -mp: if -m >2000 -mp on
//	final SettingsModelBoolean ERROR_FREE = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_ERROR_FREE, false); //if true, use -ef
//	final SettingsModelBoolean SEPERATE_PROFILE = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_SEPERATE_PROFILE, false); //if true, use -sp
//	final SettingsModelBoolean QUIET = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_QUIET, true); //if false, do not use -q

	final SettingsModelString use_FILE = new SettingsModelString(ArtNodeModel.CFGKEY_ID_PATH,"");
    final SettingsModelBoolean use_NO_MASK = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_NO_MASK, false); //if true, use -nf 0
	final SettingsModelBoolean use_MATE_PAIR = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_MATE_PAIR, false); //if true, use -mp: if -m >2000 -mp on
	final SettingsModelBoolean use_ERROR_FREE = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_ERROR_FREE, false); //if true, use -ef
	final SettingsModelBoolean use_SEPERATE_PROFILE = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_SEPERATE_PROFILE, false); //if true, use -sp
	final SettingsModelBoolean use_QUIET = new SettingsModelBoolean(ArtNodeModel.CFGKEY_USE_QUIET, false); //if false, do not use -q

    protected ArtNodeDialog() {
        super();
        
    	setHorizontalPlacement(false);
    	
    	String s = "";
    	s = (ArtNodeModel.optionalPort)? "enabled": "diasbled";
    	addDialogComponent(new DialogComponentFileChooser(use_FILE, "Testing", 0, true));

//    	addDialogComponent(new DialogComponentFileChooser(use_FILE,"Testing", 0, ".txt"));
//    	else
//    		addDialogComponent(new DialogComponentBoolean(new SettingsModelBoolean("optionalport", true), "Optional port enabled"));
    	


      	addDialogComponent(new DialogComponentNumber(length,"Read length(bp):", /*step*/ 10, /*componentwidth*/ 5));
      	addDialogComponent(new DialogComponentNumber(mean_size,"Mean size of DNA fragment:", /*step*/ 50, /*componentwidth*/ 5));     
      	addDialogComponent(new DialogComponentNumber(fold,"Fold of read coverage:", /*step*/ 1, /*componentwidth*/ 5));
      	addDialogComponent(new DialogComponentNumber(sd,"Standard deviation:", /*step*/ 1, /*componentwidth*/ 5));
      	
      	addDialogComponent(new DialogComponentBoolean(use_NO_MASK, "Turn off the masking of N region"));
      	addDialogComponent(new DialogComponentBoolean(use_MATE_PAIR, "Turn on mate pair condition"));
      	addDialogComponent(new DialogComponentBoolean(use_ERROR_FREE, "Generate the zero sequencing errors SAM file as well the regular one"));
      	addDialogComponent(new DialogComponentBoolean(use_SEPERATE_PROFILE, "Use separate quality profiles for different bases (ATGC)"
      	/*+"\n"+"NOTE: art will automatically switch to a mate-pair simulation if the given mean fragment size >= 2000"*/));
      	addDialogComponent(new DialogComponentBoolean(use_QUIET, "Turn off end of run summary"));


      	/**
         * checkboxes 
         */
//  		mutation.setEnabled(false);
//  		recombination.setEnabled(false);
//  		seed.setEnabled(false);

 
//      	if(ArtNodeModel.optionalPort) {
//      		use_FILE.addChangeListener(new ChangeListener() {
//      			public void stateChanged(ChangeEvent e) {
//      				use_FILE.setEnabled(false); 				
//      			}
//          	});
//      	}
//      	
      	use_NO_MASK.addChangeListener(new ChangeListener() {
  			public void stateChanged(ChangeEvent e) {
  				use_NO_MASK.setEnabled(true); 				
  			}
      	});
      	use_MATE_PAIR.addChangeListener(new ChangeListener() {
  			public void stateChanged(ChangeEvent e) {
  				use_MATE_PAIR.setEnabled(true); 				
  			}
      	});
      	use_ERROR_FREE.addChangeListener(new ChangeListener() {
  			public void stateChanged(ChangeEvent e) {
  				use_ERROR_FREE.setEnabled(true); 				
  			}
      	});
      	use_SEPERATE_PROFILE.addChangeListener(new ChangeListener() {
  			public void stateChanged(ChangeEvent e) {
  				use_SEPERATE_PROFILE.setEnabled(true); 				
  			}
      	});
      	use_QUIET.addChangeListener(new ChangeListener() {
  			public void stateChanged(ChangeEvent e) {
  				use_QUIET.setEnabled(true); 				
  			}
      	});
      	
      	

                    
    }
}

