list := 10 100 1000 10000

all: one_exp loop_exp

one_exp:
	./uebung1 -f ../data/oneBigexperiment.txt -p 10000 -x 1 -s 1 > ../data/one_exp_estimator.txt
	./uebung1 -f ../data/split_100_exp.txt -p 100 -x 100 -s 0 > ../data/split_100_exp_estimator.txt 
	./uebung1 -f ../data/zillion_little_exp.txt -p 1 -x 10000 -s 0 > ../data/zillion_little_exp_estimator.txt

loop_exp: 
	$(foreach i,$(list), \
		$(foreach j,$(list), \
		./uebung1 -p $(i) -x $(j) -s False > ../data/combination_pairs_$(i)_exp_$(j).txt ; \
		)\
	)

clean:
	rm -f uebung1 ../data/*