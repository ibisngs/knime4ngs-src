package de.helmholtz_muenchen.ibis.ngs.gatkraw;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GATKRaw" Node.
 * GATK node without given command walker
 *
 * @author Tim Jeske
 */
public class GATKRawNodeFactory 
        extends NodeFactory<GATKRawNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKRawNodeModel createNodeModel() {
        return new GATKRawNodeModel(0,1);
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
    public NodeView<GATKRawNodeModel> createNodeView(final int viewIndex,
            final GATKRawNodeModel nodeModel) {
        return new GATKRawNodeView(nodeModel);
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
        return new GATKRawNodeDialog();
    }

}

