package de.helmholtz_muenchen.ibis.ngs.fastqSplitter;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "FastqSplitter" Node.
 * 
 *
 * @author Michael Kluge
 */
public class FastqSplitterNodeFactory 
        extends NodeFactory<FastqSplitterNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FastqSplitterNodeModel createNodeModel() {
        return new FastqSplitterNodeModel();
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
    public NodeView<FastqSplitterNodeModel> createNodeView(final int viewIndex,
            final FastqSplitterNodeModel nodeModel) {
        return new FastqSplitterNodeView(nodeModel);
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
        return new FastqSplitterNodeDialog();
    }

}

