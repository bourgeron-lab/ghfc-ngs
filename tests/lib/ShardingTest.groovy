// Sharding decides where every published file lives; a change here moves the whole data tree.

Check.that('shards are the last two characters, last first, upper-cased') {
    assert Sharding.getShards('fam1a') == ['A', '1']
}
Check.that('separators are ignored when sharding') {
    assert Sharding.getShards('C000F2W') == Sharding.getShards('C000-F2.W_')
    assert Sharding.getShards('FAM-42') == ['2', '4']
}
Check.that('sample and family directories') {
    assert Sharding.getSampleDir('/d', 'C000F2W') == '/d/samples/W/2/C000F2W'
    assert Sharding.getFamilyDir('/d', 'FAM001') == '/d/families/1/0/FAM001'
}

Check.done()
