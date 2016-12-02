import Queue

"""
Functions to deal with taking data from queue

get_queue_item returns first item from queue
(queue default is FIFO)

get_all_queue_items returns all items in the
queue
"""

def get_queue_item(queue_object):
    try:
        item = queue_object.get(False)
    except Queue.Empty:
        return None
    return item

def get_all_queue_items(queue_object):
    pass
