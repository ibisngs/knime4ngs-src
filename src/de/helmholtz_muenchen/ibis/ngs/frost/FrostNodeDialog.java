package de.helmholtz_muenchen.ibis.ngs.frost;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentBoolean;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.DialogComponentNumber;
import org.knime.core.node.defaultnodesettings.DialogComponentStringSelection;
import org.knime.core.node.defaultnodesettings.SettingsModelBoolean;
import org.knime.core.node.defaultnodesettings.SettingsModelDoubleBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelIntegerBounded;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

/**
 * <code>NodeDialog</code> for the "Frost" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Syeda Tanzeem Haque
 */
public class FrostNodeDialog extends DefaultNodeSettingsPane {
	
	final SettingsModelString parameter = new SettingsModelString(FrostNodeModel.CFGKEY_VARY,"");
	final SettingsModelDoubleBounded mutation = new SettingsModelDoubleBounded(FrostNodeModel.CFGKEY_MUT_RATE,FrostNodeModel.DEFAULT_MUTATION_RATE,1.0, 3.0);
    final SettingsModelIntegerBounded recombination = new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_RECOMB_NUM,FrostNodeModel.DEFAULT_RECOMNATION,0, Integer.MAX_VALUE);
    final SettingsModelIntegerBounded generation = new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_GENERATION,FrostNodeModel.DEFAULT_GENERATION,0, Integer.MAX_VALUE);
    final SettingsModelIntegerBounded seed = new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_SEED,FrostNodeModel.DEFAULT_SEED,Integer.MIN_VALUE, Integer.MAX_VALUE);

	final SettingsModelBoolean use_MUT = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_MUT_RATE, false);
	final SettingsModelBoolean use_REC = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_RECOMB_NUM, false);
	final SettingsModelBoolean use_GEN = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_GENERATION, false);
	final SettingsModelBoolean use_SEED = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_SEED, false);

    /**
     * New pane for configuring the Frost node.
     */
    protected FrostNodeDialog() {
    	super();
    	

    	
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(FrostNodeModel.CFGKEY_FASTA,""), 
    			"Testing", 0, ".fa", ".fasta", ".FASTA"));
    	
    	createNewGroup("");
    	addDialogComponent(new DialogComponentStringSelection(parameter,"Select a parameter whose positions will vary at each run","Mutation","Crossover", "Denovo"));
      	
    	createNewGroup("");
    	addDialogComponent(new DialogComponentBoolean(use_MUT, "User specified mutation rate"));
      	addDialogComponent(new DialogComponentNumber(mutation,"Mutation Rate (*e-08):", /*step*/ 0.01, /*componentwidth*/ 5));
      	
    	createNewGroup("");
      	addDialogComponent(new DialogComponentBoolean(use_REC, "User specified number of crossovers"));
      	addDialogComponent(new DialogComponentNumber(recombination,"#Crossover-points:", /*step*/ 1, /*componentwidth*/ 5));
      
    	createNewGroup("");
    	addDialogComponent(new DialogComponentBoolean(use_GEN, "User specified number of generations since first Homo sapiens"));
      	addDialogComponent(new DialogComponentNumber(generation,"#Generation:", /*step*/ 500, /*componentwidth*/ 5));
      	
    	createNewGroup("");
      	addDialogComponent(new DialogComponentBoolean(use_SEED, "User specified seed to produce same positions for the other two fixed parameters"));
      	addDialogComponent(new DialogComponentNumber(seed,"Random seed:", /*step*/ 1, /*componentwidth*/ 5));
      	
//      	addDialogComponent(new DialogComponentNumber(new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_RECOMB_NUM,
//        		FrostNodeModel.DEFAULT_RECOMNATION,0, Integer.MAX_VALUE),"Recombination:", /*step*/ 1, /*componentwidth*/ 10));
//        addDialogComponent(new DialogComponentNumber(new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_SEED,
//        		FrostNodeModel.DEFAULT_SEED,Integer.MIN_VALUE, Integer.MAX_VALUE),"Seed:", /*step*/ 1, /*componentwidth*/ 15));


      /**
       * checkboxes 
       */
//		mutation.setEnabled(false);
//		recombination.setEnabled(false);
//		seed.setEnabled(false);

    	use_MUT.addChangeListener(new ChangeListener() {

			public void stateChanged(ChangeEvent e) {

				if(use_MUT.getBooleanValue()){
					mutation.setEnabled(true);
										
				}
				else{				
					mutation.setEnabled(false);
					
				}
			}
    	});
    	use_REC.addChangeListener(new ChangeListener() {

			public void stateChanged(ChangeEvent e) {

				if(use_REC.getBooleanValue()){
					recombination.setEnabled(true);
										
				}
				else{				
					recombination.setEnabled(false);
					
				}
			}
    	});
    	use_GEN.addChangeListener(new ChangeListener() {

			public void stateChanged(ChangeEvent e) {

				if(use_GEN.getBooleanValue()){
					generation.setEnabled(true);
										
				}
				else{				
					generation.setEnabled(false);
					
				}
			}
    	});
    	use_SEED.addChangeListener(new ChangeListener() {

			public void stateChanged(ChangeEvent e) {

				if(use_SEED.getBooleanValue()){
					seed.setEnabled(true);
										
				}
				else{				
					seed.setEnabled(false);
					
				}
			}
    	});
    	
    	setHorizontalPlacement(true);

    }
}

