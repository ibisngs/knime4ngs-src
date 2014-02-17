package de.helmholtz_muenchen.ibis.ngs.fastaSelector;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.FileSelector.FileSelectorNodeModel;


/**
 * <code>FastaSelectorNodeModel</code> can be used to select multiple fasta files.
 *
 * @author Michael Kluge
 */
public class FastaSelectorNodeModel extends FileSelectorNodeModel {

	public static final String OUTPUT_NAME_FASTA_FILES 	= "FastaFiles";
	public static final String FILTERTYPE_NAME 			= "FASTA";
	public static final String FILENAME_END_REGEX 		= "\\.(fa|fasta)$";
	
	@Override
	public String getNameOfOutputCol() {
		return OUTPUT_NAME_FASTA_FILES;
	}
	
	@Override
	protected String getFiletypeName() {
		return FILTERTYPE_NAME;
	}
}

