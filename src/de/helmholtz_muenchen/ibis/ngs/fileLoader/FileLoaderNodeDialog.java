package de.helmholtz_muenchen.ibis.ngs.fileLoader;

import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.knime.core.node.defaultnodesettings.DefaultNodeSettingsPane;
import org.knime.core.node.defaultnodesettings.DialogComponentFileChooser;
import org.knime.core.node.defaultnodesettings.SettingsModelString;

import de.helmholtz_muenchen.ibis.utils.IO;


/**
 * <code>NodeDialog</code> for the "FileLoader" Node.
 * 
 *
 * This node dialog derives from {@link DefaultNodeSettingsPane} which allows
 * creation of a simple dialog with standard components. If you need a more 
 * complex dialog please derive directly from 
 * {@link org.knime.core.node.NodeDialogPane}.
 * 
 * @author Tim Jeske
 */
public class FileLoaderNodeDialog extends DefaultNodeSettingsPane {

	private final SettingsModelString m_infile1 = new SettingsModelString(FileLoaderNodeModel.CFGKEY_INFILE1,"");
	private final SettingsModelString m_infile2 = new SettingsModelString(FileLoaderNodeModel.CFGKEY_INFILE2,"");
    /**
     * New pane for configuring the FileLoader node.
     */
    protected FileLoaderNodeDialog() {
    	createNewGroup("Input file (BAM/SAM, fastQ, VCF)");
    	addDialogComponent(new DialogComponentFileChooser(m_infile1, "his0_infile1", 0, ".bam|.sam","fastq|.fastq.gz|.fq|.fq.gz",".vcf|.vcf.gz|.gvcf",".csv|.tsv|.list"));

    	createNewGroup("Second optional fastQ file");
    	addDialogComponent(new DialogComponentFileChooser(m_infile2, "his0_infile2", 0, "fastq|.fastq.gz|.fq"));
    
    	m_infile1.addChangeListener(new ChangeListener() {

			@Override
			public void stateChanged(ChangeEvent e) {
				String file1 = m_infile1.getStringValue();
				file1 = IO.removeZipExtension(file1);
				if(file1.endsWith(".fastq") || file1.endsWith(".fq")) {
					m_infile2.setEnabled(true);
				} else {
					m_infile2.setEnabled(false);
					m_infile2.setStringValue("");
				}
			}
    	});
    
    }
    
}

