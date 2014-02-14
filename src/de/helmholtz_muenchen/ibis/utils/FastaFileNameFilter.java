package de.helmholtz_muenchen.ibis.utils;

import java.io.File;
import java.io.FilenameFilter;

/**
 * FilenameFilter for fasta files which accepts files with ending .fa or .fasta
 * @author Michael Kluge
 *
 */
public class FastaFileNameFilter implements FilenameFilter
{
	private static final String FASTA = ".fasta";
	private static final String FA = ".fa";

	@Override
	public boolean accept(File dir, String name) {
		return name.endsWith(FASTA) || name.endsWith(FA);
	}
}