aws s3 cp  s3://DeepMEI/regions_bam/$1 . --endpoint-url http://10.91.230.2:9020
docker exec -i -w /DeepMEI/DeepMEI_model deepAlu  bash  model_test_batch.sh -i /DeepMEI/$1
rm ${1}*
