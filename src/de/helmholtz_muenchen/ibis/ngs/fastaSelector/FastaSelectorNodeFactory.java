package de.helmholtz_muenchen.ibis.ngs.fastaSelector;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "FastaSelector" Node.
 * This Node can be used to select multiple fasta files.
 *
 * @author Michael Kluge
 */
public class FastaSelectorNodeFactory 
        extends NodeFactory<FastaSelectorNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FastaSelectorNodeModel createNodeModel() {
        return new FastaSelectorNodeModel();
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
    public NodeView<FastaSelectorNodeModel> createNodeView(final int viewIndex, final FastaSelectorNodeModel nodeModel) {
        return new FastaSelectorNodeView(nodeModel);
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
        return new FastaSelectorNodeDialog();
    }
}

