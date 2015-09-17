package de.helmholtz_muenchen.ibis.ngs.frost;

import javax.swing.JFileChooser;
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
    /**
     * final SettingsModelIntegerBounded recombination = new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_RECOMB_NUM,FrostNodeModel.DEFAULT_RECOMNATION,0, Integer.MAX_VALUE);
     */
    /**
     * final SettingsModelIntegerBounded generation = new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_GENERATION,FrostNodeModel.DEFAULT_GENERATION,0, Integer.MAX_VALUE);
     */
    final SettingsModelIntegerBounded seed = new SettingsModelIntegerBounded(FrostNodeModel.CFGKEY_SEED,FrostNodeModel.DEFAULT_SEED,Integer.MIN_VALUE, Integer.MAX_VALUE);
	final SettingsModelString bedFile = new SettingsModelString(FrostNodeModel.CFGKEY_BED_FILE, "");

	final SettingsModelBoolean use_MUT = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_MUT_RATE, false);
	/**
	 * final SettingsModelBoolean use_REC = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_RECOMB_NUM, false);
	 */
	/**
	 * final SettingsModelBoolean use_GEN = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_GENERATION, false);
	 */
	final SettingsModelBoolean use_SEED = new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_SEED, false);
    final SettingsModelBoolean use_BED_FILE	= new SettingsModelBoolean(FrostNodeModel.CFGKEY_USE_BED_FILE, false);

    /**
     * New pane for configuring the Frost node.
     */
    protected FrostNodeDialog() {
    	super();
    	

    	
    	addDialogComponent(new DialogComponentFileChooser(new SettingsModelString(FrostNodeModel.CFGKEY_FASTA,""), 
    			"Testing", 0, ".fa", ".fasta", ".FASTA"));
    	
    	createNewGroup("Parameters");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentStringSelection(parameter,"Select one to vary the positions at each run","Mutation","Crossover", "Denovo"));
    	setHorizontalPlacement(false);

    	createNewGroup("Exome simulation");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(use_BED_FILE, "Use .bed file"));
    	addDialogComponent(new DialogComponentFileChooser(bedFile, "BED_FILE", JFileChooser.OPEN_DIALOG, false, ".bed"));
		setHorizontalPlacement(false);

    	createNewGroup("User specifications");
    	setHorizontalPlacement(true);
    	addDialogComponent(new DialogComponentBoolean(use_MUT, "Mutation Rate (*e-08):"/*"Mutation rate"*/));
      	addDialogComponent(new DialogComponentNumber(mutation,"", /*step*/ 0.01, /*componentwidth*/ 5));
		setHorizontalPlacement(false);
		/**
		 * 
		 */
//
////    	createNewGroup("Recombination");
//    	setHorizontalPlacement(true);
//      	addDialogComponent(new DialogComponentBoolean(use_REC, "#Crossover-points:"/*"Number of crossovers"*/));
//      	addDialogComponent(new DialogComponentNumber(recombination,"", /*step*/ 1, /*componentwidth*/ 5));
//		setHorizontalPlacement(false);
//		
////    	createNewGroup("Generation");
//    	setHorizontalPlacement(true);
//    	addDialogComponent(new DialogComponentBoolean(use_GEN, "#Generation:"/*"Generations since first Homo sapiens"*/));
//      	addDialogComponent(new DialogComponentNumber(generation,"", /*step*/ 500, /*componentwidth*/ 5));
//		setHorizontalPlacement(false);
		/**
		 * 
		 */
//    	createNewGroup("Seed");
    	setHorizontalPlacement(true);
      	addDialogComponent(new DialogComponentBoolean(use_SEED, "Random seed:"/*"Seed to produce same positions for the other two fixed parameters"*/));
      	addDialogComponent(new DialogComponentNumber(seed,"", /*step*/ 1, /*componentwidth*/ 5));
		setHorizontalPlacement(false);

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
				mutation.setEnabled(use_MUT.getBooleanValue());
			}
    	});
    	/**
    	use_REC.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				recombination.setEnabled(use_REC.getBooleanValue());
			}
    	});
    	use_GEN.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				generation.setEnabled(use_GEN.getBooleanValue());
			}
    	});**/
    	use_SEED.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
				seed.setEnabled(use_SEED.getBooleanValue());
			}
    	});
    	use_BED_FILE.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent e) {
					bedFile.setEnabled(use_BED_FILE.getBooleanValue());
				}
		});
    	
//    	setHorizontalPlacement(true);

    }
}

