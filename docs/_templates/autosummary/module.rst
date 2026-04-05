{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}

{% if modules %}
Submodules
----------

.. autosummary::

{% for item in modules %}
   {{ item }}
{%- endfor %}
{% endif %}

{% if attributes %}
Module Attributes
-----------------

.. autosummary::

{% for item in attributes %}
   {{ item }}
{%- endfor %}
{% endif %}

{% if functions %}
Functions
---------

.. autosummary::

{% for item in functions %}
   {{ item }}
{%- endfor %}

Function Reference
------------------

{% for item in functions %}
.. autofunction:: {{ fullname }}.{{ item }}

{% endfor %}
{% endif %}

{% if classes %}
Classes
-------

.. autosummary::

{% for item in classes %}
   {{ item }}
{%- endfor %}

Class Reference
---------------

{% for item in classes %}
.. autoclass:: {{ fullname }}.{{ item }}
   :members:
   :show-inheritance:

{% endfor %}
{% endif %}

{% if exceptions %}
Exceptions
----------

.. autosummary::

{% for item in exceptions %}
   {{ item }}
{%- endfor %}

Exception Reference
-------------------

{% for item in exceptions %}
.. autoexception:: {{ fullname }}.{{ item }}

{% endfor %}
{% endif %}
