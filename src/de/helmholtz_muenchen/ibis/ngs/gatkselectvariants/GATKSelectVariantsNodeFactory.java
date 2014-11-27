package de.helmholtz_muenchen.ibis.ngs.gatkselectvariants;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GATKSelectVariants" Node.
 * 
 *
 * @author 
 */
public class GATKSelectVariantsNodeFactory 
        extends NodeFactory<GATKSelectVariantsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKSelectVariantsNodeModel createNodeModel() {
        return new GATKSelectVariantsNodeModel(1,1);
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
    public NodeView<GATKSelectVariantsNodeModel> createNodeView(final int viewIndex,
            final GATKSelectVariantsNodeModel nodeModel) {
        return new GATKSelectVariantsNodeView(nodeModel);
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
        return new GATKSelectVariantsNodeDialog();
    }

}

