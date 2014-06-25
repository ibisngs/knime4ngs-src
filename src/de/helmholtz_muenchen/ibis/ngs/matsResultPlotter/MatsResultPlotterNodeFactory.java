package de.helmholtz_muenchen.ibis.ngs.matsResultPlotter;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "MatsResultPlotter" Node.
 * 
 *
 * @author Michael Kluge
 */
public class MatsResultPlotterNodeFactory 
        extends NodeFactory<MatsResultPlotterNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public MatsResultPlotterNodeModel createNodeModel() {
        return new MatsResultPlotterNodeModel();
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
    public NodeView<MatsResultPlotterNodeModel> createNodeView(final int viewIndex,
            final MatsResultPlotterNodeModel nodeModel) {
        return new MatsResultPlotterNodeView(nodeModel);
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
        return new MatsResultPlotterNodeDialog();
    }

}

