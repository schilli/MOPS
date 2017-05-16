try:
    # python 3 ?
    from MDorder.MDorder import save_corr, load_corr, bondvec_corr, bondvec_corr_batch
    from MDorder.OrderParameter import OrderParameter
    #from MDorder.corrFunction import corrFunction
except ImportError:
    # python 2 ?
    from MDorder import save_corr, load_corr, bondvec_corr, bondvec_corr_batch
    from OrderParameter import OrderParameter 
    #from corrFunction import corrFunction

__all__ = []
