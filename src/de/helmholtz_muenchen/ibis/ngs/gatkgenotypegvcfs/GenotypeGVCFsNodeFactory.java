package de.helmholtz_muenchen.ibis.ngs.gatkgenotypegvcfs;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GenotypeGVCFs" Node.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class GenotypeGVCFsNodeFactory 
        extends NodeFactory<GenotypeGVCFsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GenotypeGVCFsNodeModel createNodeModel() {
        return new GenotypeGVCFsNodeModel(1,1);
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
    public NodeView<GenotypeGVCFsNodeModel> createNodeView(final int viewIndex,
            final GenotypeGVCFsNodeModel nodeModel) {
        return new GenotypeGVCFsNodeView(nodeModel);
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
        return new GenotypeGVCFsNodeDialog();
    }

}

