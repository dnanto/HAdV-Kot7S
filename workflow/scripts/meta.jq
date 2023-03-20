{
    "species": input_filename | split("/")[2], 
    "accession": .accession, 
    "sourceDatabase": .sourceDatabase, 
    "virus.taxId": .virus.taxId, 
    "virus.organismName": .virus.organismName, 
    "isolate.name": .isolate.name, 
    "isolate.collectionDate": .isolate.collectionDate, 
    "location.geographicLocation": .location.geographicLocation, 
    "location.geographicRegion": .location.geographicRegion, 
    "host.taxId": .host.taxId
}
