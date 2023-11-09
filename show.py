import networkx as nx
import matplotlib.pyplot as plt


def main():

    # read stdin while not EOF
    while True:
        try:
            line = input().split()

            if len(line) == 0:
                continue

            if line[0][0] != "#":
                continue

            assert len(line) >= 4

            graph_name, ns, ms, ls = line[:4]
            assert ls == "0"

            line = line[4:]


            G = nx.Graph()
            
            for i in range(0, int(ns)):
                G.add_node(i)

            for i in range(0, int(ms)):
                a, b = line[:2]
                G.add_edge(int(a), int(b))
                line = line[2:]
                
            
            # # show the graphs


            plt.figure()
            plt.title(graph_name)
            nx.draw(G, with_labels=True)
            plt.show()



        # if EOF, break
        except EOFError:
            break




if __name__ == "__main__":
    main()