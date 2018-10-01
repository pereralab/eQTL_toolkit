import gzip 

def threshold(f):
        eqtl = {}
        threshold = open(f)
        out = open(f + ".significant.txt", "w")
        for line in threshold:
                items = line.split()
                if float(items[2]) <= 0.05:
                        out.write(line)
                        eqtl[items[0]] = float(items[1])
        out.close()                        
        return eqtl                


def call(egene, nominal, out):
        eqtl_f = open(nominal)
        eqtl_out = open(out,"w")
        for line in eqtl_f:
                items = line.split()
                
                if float(items[3]) <= egene[items[0]]:
                        eqtl_out.write(line)
        eqtl_out.close()  
