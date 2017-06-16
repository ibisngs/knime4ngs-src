package de.helmholtz_muenchen.ibis.ngs.trimgalore;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.HTExecutorNode.HTENodeView;

/**
 * <code>NodeFactory</code> for the "TrimGalore" Node.
 * Trim Galore! is a wrapper script to automate quality and adapter trimming as well as quality control, with some added functionality to remove biased methylation positions for RRBS sequence files (for directional, non-directional (or paired-end) sequencing).
 *
 * @author Paul Hager
 */
public class TrimGaloreNodeFactory 
        extends NodeFactory<TrimGaloreNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public TrimGaloreNodeModel createNodeModel() {
        return new TrimGaloreNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<TrimGaloreNodeModel> createNodeView(final int viewIndex,
            final TrimGaloreNodeModel nodeModel) {
        return new HTENodeView<TrimGaloreNodeModel>(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new TrimGaloreNodeDialog();
    }

}

