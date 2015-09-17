package de.helmholtz_muenchen.ibis.ngs.limma;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeFactory</code> for the "Limma" Node.
 * 
 *
 * @author Michael Kluge
 */
public class LimmaNodeFactory 
        extends NodeFactory<LimmaNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public LimmaNodeModel createNodeModel() {
        return new LimmaNodeModel();
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
    public NodeView<LimmaNodeModel> createNodeView(final int viewIndex,
            final LimmaNodeModel nodeModel) {
        return new RNodeView<LimmaNodeModel>(nodeModel);
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
        return new LimmaNodeDialog();
    }

}

