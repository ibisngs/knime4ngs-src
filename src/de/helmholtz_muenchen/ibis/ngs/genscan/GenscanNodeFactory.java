package de.helmholtz_muenchen.ibis.ngs.genscan;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Genscan" Node.
 * 
 *
 * @author Jan
 */
public class GenscanNodeFactory 
        extends NodeFactory<GenscanNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GenscanNodeModel createNodeModel() {
        return new GenscanNodeModel();
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
    public NodeView<GenscanNodeModel> createNodeView(final int viewIndex,
            final GenscanNodeModel nodeModel) {
        return new GenscanNodeView(nodeModel);
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
        return new GenscanNodeDialog();
    }

}

