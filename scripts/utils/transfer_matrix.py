class TransferMatrix(object):
    """
    Transfer Matrix uses 
    """
    def __init__(self, M):
      
     @classmethod
     def calculate_m(cls, x_in, x_out):
        """
        Calculate a  transfer matrix m based on a number of points x_in and x_out
        """
        if x_in.shape != x_out.shape:
            raise ValueError("Need equal shape of input and output points")
        if x_in.shape[0] != 2*x_in.shape[1]:
            raise ValueError("Should have twice as many points as dimensions")
        dx = 
        for i in range(x_in.shape[0])
        

    def print()